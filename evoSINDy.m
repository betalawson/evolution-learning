function K = evoSINDy(X_data, X_dash_data, F, options)
% 
% K = evoSINDy(X_data, X_dash_data, F, options)
%
% This function takes a series of observations regarding feature frequency,
% except expressed in terms of how the current state of the system relates
% to its rate of change, X vs X_dash, and determines coefficients for a
% library of functions, F, that are though to in some combination specify
% this relationship. Owing to learning of a fitness function rather than a
% right hand side function, to access the SINDy framework we must also
% transform the fitness functions F into their replicator-transformed
% equivalents. That is, we learn
%    X' = X . ( sum[k_i F_i(X)] - sum[k_i F_i(X)^T X] e ),
% where 'e' is the vector of all ones and each element F_i(X) is one of the
% functions in the full provided set of functions 'F'. 
% 
% The replicator transform converts this into the typical SINDy learning 
% problem:
%    X' = sum[ k_i G_i(X) ].
% 
% The user can optionally supply options regarding whether to force 
% positivity of learned coefficients, and the level of shrinkage to the 
% null model that is applied. If these options are not supplied, default 
% options are used.
%
% -Inputs-
%
%         X_data: An N_s by N_obs matrix containing N_obs observations of
%                 the system state, as defined by the values of its N_s
%                 species
%
%    X_dash_data: An N_s by N_obs matrix containing the observations of the
%                 time derivative (or estimations of such) of the system's
%                 species, each corresponding to the observations of system
%                 state in input 'X_data'
%
%              F: A cell array specifying the library functions, each
%                 element of the array a function of form f' = @(x) that
%                 takes as input the system state as an N_s by 1 vector and
%                 outputs N_s by 1 contributions to the derivative
%
%        options: Used here to specify whether to enforce positivity of
%                 constants, and also to provide the level of shrinkage
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INITIAL PREPARATION

% Read out the number of species and number of observations
[N_spec, N_obs] = size(X_data);
N_funs = length(F);

% If options not provided, get defaults
if nargin < 4
    options = [];
end
options = addDefaultOptions(options);

% Transform the set of input functions into type I replicator functions
F = convertToReplicatorLibrary(F);


%%% CONSTRUCTION OF SINDy MATRIX SYSTEM

% Create the "A matrix":
% Library functions go across columns, and rows correspond to the
% observations in the data (with the different dependent species stacked on
% top of one another)

% Initialise storage
XI = zeros(N_spec*N_obs,N_funs);

% Loop over each library function and evaluate it for the system states in 
% the data to build the matrix
for k = 1:N_funs
    
    % Evaluate the current function at all the X values in the data
    % (Each row is a dependent species, each column is a datapoint)
    F_eval = F{k}(X_data);
    
    % Convert these function evaluations into a single column to go in the
    % matrix for this reaction (library function)
    XI(:,k) = F_eval(:);

end

% Create the "b vector" using the derivative data
b = X_dash_data(:);


%%% SOLVE THE MATRIX SYSTEM

% If using shrinkage, add "pseudo-observations" that penalise coefficients
if options.shrinkage > 0
    XI = [XI; sqrt(options.shrinkage) * eye(N_funs)];
    b = [b; zeros(N_funs,1)]; 
end
    
% Use minimum norm solution if unconstrained, otherwise use built-in
% MATLAB function for positive least squares solution
if options.force_positive
    K = lsqnonneg(XI' * XI,XI' * b);
else
    K = lsqminnorm(XI' * XI,XI' * b);
end