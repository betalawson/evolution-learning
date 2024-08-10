function K = SINDy(X_data, X_dash_data, F, options)
% 
% K = SINDy(X_data, X_dash_data, F, options)
%
% This function takes a series of observations regarding how the state of a
% system relates to its rate of change, as specified by inputs X_data and
% X_dash_data, and attempts to determine the coefficients for a library of
% functions, F, that specify this relationship (and hence specify the set
% of ordinary differential equations that govern that system). The user can
% also optionally specify that coefficients must be positive.
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
% Optional arguments:
%
%     interactions - a matrix of true/false values, indicating which
%     species are allowed to interact with which other species in
%     constructing the list of reactions
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
options = imbueELDefaults(options);


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