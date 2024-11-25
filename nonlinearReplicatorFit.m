function K = nonlinearReplicatorFit(t_data, X_data, library, selection_type, options)
%
%     K = nonlinearReplicatorFit(X_data, X_dash_data, F, selection_type)
%     K = nonlinearReplicatorFit(X_data, X_dash_data, F, selection_type, options)
%
% This function takes a series of observations of feature frequency data
% over time, as specified by inputs t_data and X_data, and attempts to
% determine the coefficients for a library of functions, F, that when taken
% together to determine the fitness, produce replicator trajectories that
% match that data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set all options that weren't pre-specified to their defaults
if nargin < 4
    options = [];
end
options = addDefaultOptions(options);

% Use options to decide whether to supress common nlinfit warnings
if options.nlin_silent
    warning('off','all');
end 

% Use the built-in nonlinear fitter,    BETA = nlinfit(X,Y,MODELFUN,BETA0)
% Covariates (X) are the timepoints where data was recorded, t_data.
% Responses (Y) are the observed system states at those timepoints.
% MODELFUN is defined below, and simulates the replicator equation.
% Initial condition (BETA0) is taken as all ones.

% Set "base" model according to the type of selection - this is mostly just
% to match the approaches taken for gradient matching
N_funs = length(library.F);
if selection_type == 2
    K0 = 1/sqrt(N_funs) * ones(N_funs,1);
else
    K0 = zeros(N_funs,1);
end

% Append a value of zero as "data" that is used to punish parameters
% differing from the base model
Yfit = [X_data(:); 0];

% Grab out the initial point in the data as the starting point
X0 = X_data(:,1);

% If forcing positive, optimise by nonlinear fitting in square root space
if options.force_positive
        
    % Create fitting function with square transform to force positive and
    % append "cost" of large values for coefficients
    fitfun = @(sqrtK,t) [ evaluateReplicator(t,sqrtK.^2,X0,library,selection_type); sqrt(options.shrinkage)*norm(sqrtK.^2 - K0) ];
    
    % Fit the model
    sqrtK = nlinfit( t_data, Yfit, fitfun, ones(N_funs,1));
    K = sqrtK.^2;
    
else
    
    % Create fitting function with square transform to force positive and
    % append "cost" of large values for coefficients
    fitfun = @(K,t) [ evaluateReplicator(t,K,X0,library,selection_type); sqrt(options.shrinkage)*norm(K - K0) ];
    
    % Fit the model
    K = nlinfit( t_data, Yfit, fitfun, ones(N_funs,1));
    
end

% If warnings were turned off, re-enable them to protect user's session
if options.nlin_silent
    warning('on','all');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
function X = evaluateReplicator(t,K,X0,library,selection_type)

%%% PARAMETER DEFINITIONS

% Number of points to use for replicator trajectories
N_simpts = 10001;
% Initial tolerance to trial replicator solves
init_tol = 1e-5;
% Minimum tolerance for replicator solves before abandoning
min_tol = 1e-11;
% Tolerance decrease factor (scales tolerance upon integration failure)
tol_factor = 0.01;
% Flag for solving using either a standard or DAE formulation
solve_as_DAE = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare a mass matrix to confine results to probability simplex
N_feat = length(X0);
M = eye(N_feat);
M(N_feat,N_feat) = 0;
running = true;
success = false;

% Disable warnings as we will try many tolerances
warning('off','MATLAB:ode15s:IntegrationTolNotMet');
        
% Prepare timepoints for simulation
timepoints = linspace(0, max(t), N_simpts);

% Initialise tolerance
tol = init_tol;

% Generate the fitness function for the input functions and coefficients
%fitness = @(X) evaluateLibrary(X, F, K);
library = setfield(library,'coeffs',K);
[A,Q] = libraryToPayoff(library,N_feat);
fitness = matricesToFitness(A,Q);

% Solve the replicator ODE
while running
    
    % Attempt a DAE solve if the flag to do so is set
    DAE_failed = false;
    if solve_as_DAE
        
        % Attempt solve at current tolerance
        try
            odeoptions = odeset('Mass',M,'AbsTol',tol,'RelTol',tol);
            odesol = ode15s( @(t,X) replicatorRHS(X,fitness,selection_type,true), timepoints, X0, odeoptions );
        catch
            DAE_failed = true;
        end
        
        % If simulation did not make it to final time, mark failure
        if ~DAE_failed && abs(odesol.x(end) - max(t)) > 1e-10
            DAE_failed = true;
        end
        
    end
    
    % If DAE solve failed, or not trying, solve using standard ODE
    if ~solve_as_DAE || DAE_failed
        odeoptions = odeset('AbsTol',tol,'RelTol',tol);
        odesol = ode15s( @(t,X) replicatorRHS(X,fitness,selection_type,false), timepoints, X0, odeoptions);
    end
    
    % Terminate if integration successful
    if abs(odesol.x(end) - max(t)) <= 1e-10
        running = false;
        success = true;
    end
    
    % Terminate without success if tolerance too small
    if tol <= min_tol && ~success
        running = false;
        fprintf('\n WARNING: There was a failure to solve the ODE even with strictest tolerance! Output will be truncated.\n');
    end
    
    % Otherwise, decrease tolerance and try again
    tol = tol * tol_factor;
    
end  

% Re-enable warnings to not bother the user's session
warning('on','MATLAB:ode15s:IntegrationTolNotMet');

% If solved successfully, return results
if success
    X = deval(odesol,t);
    X = X(:);
% Otherwise, return junk data to tell nonlinear solver this is a bad spot
else
    X = 1e10 * zeros(length(X0)*length(t),1);
end

end