function [fitfun, Flib, K, texts] = evolutionLearning(data, selection_type, options)
%
%   [X_dash_data, Flib, K, texts] = evolutionLearning(data, selection_type)
%   [X_dash_data, Flib, K, texts] = evolutionLearning(data, selection_type, options)
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
%         fitfun - The functions for F(x) and F'(x) learned by gradient
%                  estimation.
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


%%% DERIVATIVE ESTIMATION

% Only perform if data wasn't supplied
if isempty(options.Xdash_data)

    % Initialise the data storage for derivative information
    Xdash_data = cell(1,N_rep);

    % Derivative estimation is performed separately for each replicate in data
    N_obs = 0;
    for k = 1:N_rep
        
        % Prepare storage
        N_obs_here = length(data{k}.t);
        Xdash_data{k} = zeros(N_feat, N_obs_here);
        % Keep a running total of number of observations for later
        N_obs = N_obs + N_obs_here;
    
        % Fit function and estimate the derivative for each feature
        f_funs = cell(N_feat,1);
        fdash_funs = cell(N_feat,1);
        for n = 1:N_feat    
            % Fit a function
            [f,fdash] = fitFunction(data{k}.t, data{k}.X(n,:),options);
            % Evaluate and store derivative estimates
            Xdash_data{k}(n,:) = fdash(data{k}.t);
            % Store the function fits themselves (as function handles)
            f_funs{n} = f;
            fdash_funs{n} = fdash;
        end    
    end

    % Create a single evaluatable function that returns a vector for f and f'
    fitfun.f_funs = f_funs;
    fitfun.fdash_funs = fdash_funs;
    
    
% Otherwise just use the provided data and return a NaN for fit function
else
    
    if ~iscell(options.Xdash_data)
        Xdash_data{1} = options.Xdash_data;
    else
        Xdash_data = options.Xdash_data;
    end
    N_obs = sum( cellfun( @(x) length(x), Xdash_data ) );
    fitfun = NaN;
end



%%% EQUATION LEARNING

% Before data compilation, prepare storage
X_data_all = zeros(N_feat, N_obs);
Xdash_data_all = zeros(N_feat, N_obs);

% Loop over each replicate and add its data to the total observation list
start_loc = 0;
for k = 1:N_rep

    % Define location in observation list for these replicates
    N_obs_here = length(data{k}.t);
    loc = start_loc + (1:N_obs_here);
    start_loc = start_loc + N_obs_here;
    % Store the data in the master datalist in this location
    X_data_all(:,loc) = data{k}.X;
    Xdash_data_all(:,loc) = Xdash_data{k};
    
end

% Regularise derivative estimates to ensure they correctly sum to zero, if
% that option was supplied
if options.evolutionary_derivative
    Xdash_data_all = Xdash_data_all - mean(Xdash_data_all);
end

% Prepare a basic library of polynomials to given order unless specifically
% targetting a symmetric payoff matrix
if ~options.symmetric_payoff
    [Flib, texts] = createPolynomialLibrary(N_feat, options.library_orders, options.names, options.interactions);
else
    [Flib, texts] = createSymmetricPayoffLibrary(N_feat, options.names, options.interactions);
end

% Type I selection allows standard SINDy after transforming library
if selection_type == 1
    
    % Convert this to a library of replicator-style functions
    Frep = convertToReplicatorLibrary(Flib);
    
    % Use standard SINDy on the master dataset
    K = SINDy(X_data_all, Xdash_data_all, Frep, options);
    
% Type II selection requires the nonlinear solver
else
    K = nonlinearReplicatorFit(X_data_all, Xdash_data_all, Flib, options);
end
    
    