function [fitfuns, Flib, K, texts] = evolutionLearning(data, selection_type, options)
%
%   [fitfuns, Flib, K, texts] = evolutionLearning(data, selection_type)
%   [fitfuns, Flib, K, texts] = evolutionLearning(data, selection_type, options)
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
%        fitfuns - The functions for F(x) and F'(x) learned by gradient
%                  estimation, one struct with fields F and Fdash per
%                  replicate
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

% Loop over the number of data replicates provided, derivative estimation
% is performed (if necessary) separately for each data replicate
N_obs = zeros(1,N_rep);
fitfuns = cell(1,N_rep);
for k = 1:N_rep
    
    % Store how many observations there are in this replicate
    N_obs(k) = size(data{k}.X,2);
    
    % If derivative data was not provided, estimate it via functional fit
    if ~isfield(data{k},'Xdash')
        
        % Use the function fitter for this set of data
        [F, Fdash] = constructSimplexFit(data{k}, options);
        
        % Evaluate and store the derivative estimates it generates
        data{k}.Xdash = Fdash(data{k}.t);
        
        % If using the smoothed function as the X data, overwrite it here
        if options.use_smoothed
            data{k}.X = F(data{k}.t);
        end
        
        % Store the fit functions here (can be useful for later plotting)
        fitfuns{k} = struct('F',F,'Fdash',Fdash);
        
    else
        
        % No functional fit was performed so just store dummy data
        fitfuns{k} = struct('F',NaN,'Fdash',NaN);
        
    end
    
end
        

%%% EQUATION LEARNING

% Before data compilation, prepare storage
X_data_all = zeros(N_feat, sum(N_obs));
Xdash_data_all = zeros(N_feat, sum(N_obs));

% Loop over each replicate and add its data to the total observation list
start_loc = 0;
for k = 1:N_rep

    % Define location in observation list for these replicates
    loc = start_loc + (1:N_obs(k));
    % Store the data in the master datalist in this location
    X_data_all(:,loc) = data{k}.X;
    Xdash_data_all(:,loc) = data{k}.Xdash;
    % Advance start location in preparation for next loop
    start_loc = start_loc + N_obs(k);
    
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
    
    