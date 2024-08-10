function K = nonlinearReplicatorFit(X_data, X_dash_data, F, options)
%
%     K = nonlinearReplicatorFit(X_data, X_dash_data, F)
%     K = nonlinearReplicatorFit(X_data, X_dash_data, F, options)
%
% This function takes a series of observations regarding how the state of a
% system relates to its rate of change, as specified by inputs X_data and
% X_dash_data, and attempts to determine the coefficients for a library of
% functions, F, that specify this relationship indirectly as contributions
% to the fitness in a replicator equation with type II selection. This
% requires a non-linear solve, as opposed to type I which is linear and
% does not require usage of this function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If options are not provided, use default options
if nargin < 4
    options = [];
end
options = imbueELDefaults(options);

% Use the built-in nonlinear fitter,    BETA = nlinfit(X,Y,MODELFUN,BETA0)
% Observations (X) are the system state observations, X_data. These are
% kept in matrix form but converted to vectors in modelfun.
% Outputs (Y) are the system rates of change, X_dash_data, as a vector.
% MODELFUN is defined below, the replicator with type II selection.
% Initial condition (BETA0) is taken as all ones.

% If forcing positive, optimise by nonlinear fitting in logspace
if options.force_positive
    
    % Create fitting function with square transform to force positive and
    % appended "cost" of large values for coefficients
    fitfun = @(K,X) [ modelfun(X,F,K.^2); options.shrinkage*norm(K.^2)^2 ];
    % Append a value of zero as "data" to punish this cost
    Yfit = [X_dash_data(:); 0];
    % Fit the model    
    sqrtK = nlinfit(X_data',Yfit,fitfun,ones(length(F),1));
    K = sqrtK.^2;
    
% Otherwise, optimise using standard nonlinear fitting
else
    
    % Create fitting function with exponential transform to force positive
    % and appended "cost" of large values for coefficients
    fitfun = @(K,X) [ modelfun(X,F,K); options.shrinkage*norm(K)^2 ];
    % Append a value of zero as "data" to punish this cost
    Yfit = [X_dash_data(:); 0];
    % Fit the model
    K = nlinfit(X_data',Yfit,fitfun,ones(length(F),1));
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subfunction that evaluates the right hand side of the type II selection
% replicator for all observations:
%       dx/dt =  x o ( f / [f^T x] - 1 )
function y_model = modelfun(X,F,k)

% Evaluate the right hand side as described
obs_RHS = X' .* ( normalised_f(X',F,k) - ones(size(X')) );
% Append data for all model species into a list for fitting purposes
y_model = obs_RHS(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subfunction that evaluates the normalised fitness values obtained given
% some current fitness function (defined by library F and coefficients k)
% and current state x. Providing a matrix X allows vectorised calculation
% across many states
function f = normalised_f(X,F,k)

% Initialise the fitness vector and normalising constant
fitvec = 0;
const = 0;

% Loop over each library function
for j = 1:length(F)
    
    % Library terms contribute directly to fitness
    fitvec = fitvec + k(j) * F{j}(X);
    % Constant is the sum of all f^T x terms
    const = const + k(j) * sum(X .* F{j}(X));
    
end

% Calculate final fitness term f / [f^T x] for output
f = fitvec ./ const;

end