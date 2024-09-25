function [Flib, K, texts] = LSevolutionLearning(data, selection_type, options)
%
%   [Flib, K, texts] = LSevolutionLearning(data, selection_type)
%   [Flib, K, texts] = LSevolutionLearning(data, selection_type, options)
%
% This function applies the evolutionary learning approach to attempt to
% identify the selective strength that is present in provided trajectory
% data. 
%
%  INPUTS:
%
%           data - The observed data, a trajectory of feature frequencies 
%                  over time. Specify data in the form of a struct with 
%                  fields "t" (vector of times) and "X" (matrix of feature
%                  frequencies at those timepoints). Or for multiple 
%                  trajectories, supply a cell array where each element 
%                  contains such a struct.
%
% selection_type - A value of "1" or "2", specifying:
%                      1: Type I selection (additive normalisation)
%                      2: Type II selection (multiplicative normalisation)
%                  Note that Type I selection is a linear solve
%
%        options - A struct of options that the user would like to specify,
%                  otherwise the default values found in imbueELDefaults.m
%                  will be used. Available options to change are listed
%                  there.
%
%
%  OUTPUTS:
%
%
%           Flib - The library of functions that was used for fitting.
%                  These are components of fitness F(x) = sum( k_i f_i(x) )
%
%              K - Vector of coefficients associated with the library Flib
%
%          texts - Names of library functions used for outputting to user
%                  or to assist in interpretation of the output.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% PRELIMINARIES

% Ensure the data is a cell array
if ~iscell(data)
    data = {data};
end

% Read out dimensions of various data pieces
N_rep = length(data);
N_feat = size(data{1}.X,1);

% Set all options that weren't pre-specified to their defaults
if nargin < 3
    options = imbueELDefaults([], N_feat);
else
    options = imbueELDefaults(options, N_feat);
end

% Use options to decide whether to supress common nlinfit warnings
if options.nlin_silent
    warning('off','all');
end 


%%% EQUATION LEARNING

% Loop over each dataset and count observations (for initialisation)
N_obs = zeros(1,N_rep);
for k = 1:N_rep
    N_obs(k) = length(data{k}.t);
end

% Now initialise data storage
t_data_all = zeros(1, sum(N_obs));
X_data_all = zeros(N_feat, sum(N_obs));

% Now loop again and store each dataset
start_loc = 0;
for k = 1:N_rep
      
    % Add this dataset's observations
    t_data_all(start_loc+(1:N_obs(k))) = data{k}.t;
    X_data_all(:,start_loc+(1:N_obs(k))) = data{k}.X;
    % Shift start location
    start_loc = start_loc + N_obs(k);
    
end

% Prepare a basic library of polynomials to given order unless specifically
% targetting a symmetric payoff matrix
if ~options.symmetric_payoff
    [Flib, texts] = createPolynomialLibrary(N_feat, options.library_orders, options.names, options.interactions);
else
    [Flib, texts] = createSymmetricPayoffLibrary(N_feat, options.names, options.interactions);
end

% Loop over each replicate and add its data to the total observation list
start_loc = 0;
for k = 1:N_rep

    % Define location in observation list for these replicates
    N_obs_here = length(data{k}.t);
    loc = start_loc + (1:N_obs_here);
    start_loc = start_loc + N_obs_here;
    % Store the data in the master datalist in this location
    X_data_all(:,loc) = data{k}.X;
    
end

% Set "base" model according to selection type
N_funs = length(Flib);
K0 = 1/sqrt(N_funs) * ones(N_funs,1);

% Append a value of zero as "data" to punish this cost
Yfit = [X_data_all(:); 0];

% If forcing positive, optimise by nonlinear fitting in square root space
if options.force_positive
        
    % Create fitting function with square transform to force positive and
    % append "cost" of large values for coefficients
    fitfun = @(sqrtK,t) [ fullRepRHS(t,sqrtK.^2,X_data_all(:,1),Flib,selection_type); sqrt(options.shrinkage)*norm(sqrtK.^2 - K0) ];
    % Fit the model
    sqrtK = nlinfit( t_data_all, Yfit, fitfun, K0 );
    K = sqrtK.^2;
    
else
    
    % Create fitting function with square transform to force positive and
    % append "cost" of large values for coefficients
    fitfun = @(K,t) [ fullRepRHS(t,K,X_data_all(:,1),Flib,selection_type); sqrt(options.shrinkage)*norm(K - K0) ];
    % Fit the model
    K = nlinfit( t_data_all, Yfit, fitfun, K0 );
    
end

% If warnings were turned off, re-enable them to protect user's session
if options.nlin_silent
    warning('on','all');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
function X = fullRepRHS(t,K,X0,Flib,selection_type)

% Prepare a mass matrix to confine results to probability simplex
N_feat = length(X0);
M = eye(N_feat);
M(N_feat,N_feat) = 0;
% Initialise ODE solve tolerance and most strict tolerance
tol = 1e-3;
min_tol = 1e-10;
running = true;
success = false;

% Disable warnings as we will try many tolerances
warning('off','MATLAB:ode15s:IntegrationTolNotMet');
        
% Solve the replicator ODE
while running
    
    % Attempt solve using the current tolerance
    odeoptions = odeset('Mass',M,'AbsTol',tol,'RelTol',tol);
    odesol = ode15s( @(t,X) replicatorRHS(X,@(Xf) libraryRHS(Xf,Flib,K),selection_type), linspace(0,max(t),10001), X0, odeoptions );
    
    % Terminate with success if integration successful
    if abs(odesol.x(end) - max(t)) < 1e-10
        running = false;
        success = true;
    end
    
    % Terminate with failure if tolerance too small
    if tol <= min_tol && ~success
        running = false;
        fprintf('\n WARNING: There was a failure to solve the ODE even with strictest tolerance! Output will be truncated.\n');
    end
        
    % Otherwise, decrease tolerance and try again
    tol = tol * 0.01;
    
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