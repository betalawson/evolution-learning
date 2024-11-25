function [library, fitfuns] = evolutionLearning(data, selection_type, options)
%
%   [library, fitfuns] = evolutionLearning(data, selection_type)
%   [library, fitfuns] = evolutionLearning(data, selection_type, options)
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
%                  otherwise the default values found in the function
%                  addDefaultOptions.m will be used. Edit that function to
%                  adjust defaults.
%
%
%  OUTPUTS:
%
%         library - A library object that represents the selective dynamics
%                   learned in response to the data. Library objects are
%                   structs with fields:
%                     F: cell array of library functions, @(x) F(x)
%                 texts: text description of each library function
%                coeffs: coefficients associated with each library function
%
%         fitfuns - The functions for F(x) and F'(x) learned by gradient
%                   estimation, one struct with fields F and Fdash per
%                   provided dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% PRELIMINARIES

% Ensure the data is a cell array
if ~iscell(data)
    data = {data};
end

% Read out dimensions of various data pieces
N_data = length(data);
N_feat = size(data{1}.X,1);

% Set all options that weren't pre-specified to their defaults
if nargin < 3
    options = addDefaultOptions([], N_feat);
else
    options = addDefaultOptions(options, N_feat);
end

% Prepare the function library to learn, currently polynomial libraries
% are supported (optionally, symmetric-constrained 1st order polynomials)
if ~options.symmetric_payoff
    library = createPolynomialLibrary(N_feat, options.library_orders, options.names, options.interactions);
else
    library = createSymmetricPayoffLibrary(N_feat, options.names, options.interactions);
end


%%% SPLIT HERE BASED ON METHOD

switch lower(options.regression_type)
    
    case {'gradient','gm','gradient-matching','gradient_matching','grad'}
        
        %%% ESTIMATE DERIVATIVES AS REQUIRED
        
        % Loop over each data replicate separately, estimating the
        % derivative if derivative data was not provided
        N_obs = zeros(1,N_data);
        fitfuns = cell(1,N_data);
        for k = 1:N_data
            
            % Store how many observations there are in this replicate
            N_obs(k) = size(data{k}.X,2);
            
            % If derivative data was not provided, estimate it via functional fit
            if ~isfield(data{k},'Xdash')
                
                % Use the function fitter for this set of data
                [F, Fdash] = fitSimplexConstrainedFunction(data{k}, options);
                
                % Evaluate and store the derivative estimates it generates
                data{k}.Xdash = Fdash(data{k}.t);
                
                % Store the fit functions here (can be useful for later plotting)
                fitfuns{k} = struct('F',F,'Fdash',Fdash);
                
            else
                
                % No functional fit was performed so just store dummy data
                fitfuns{k} = struct('F',NaN,'Fdash',NaN);
                
            end
            
        end
        
        
        %%% GATHER TOGETHER DERIVATIVE DATA
        
        % Before data compilation, prepare storage
        X_data_all = zeros(N_feat, sum(N_obs) );
        Xdash_data_all = zeros(N_feat, sum(N_obs) );

        % Loop over each data replicate to compile the total dataset
        start_loc = 0;
        for k = 1:N_data

            % Define location in observation list for these replicates
            loc = start_loc + (1:N_obs(k));
            % Store the data in the master datalist in this location
            X_data_all(:,loc) = data{k}.X;
            Xdash_data_all(:,loc) = data{k}.Xdash;
            % Advance start location in preparation for next loop
            start_loc = start_loc + N_obs(k);
    
        end

        % Numerically ensure all derivative data sums to zero, as required
        Xdash_data_all = Xdash_data_all - mean(Xdash_data_all);
        
        
        %%% USE DERIVATIVE DATA TO ESTIMATE LIBRARY COEFFICIENTS
        
        % Method used depends on the type of selection
        if selection_type == 1   
            % Type I allows "evolutionary SINDy", linear after transforming
            library.coeffs = evoSINDy(X_data_all, Xdash_data_all, library.F, options);
        else
            % Type II requires using the nonlinear fitter
            library.coeffs = nonlinearGradientFit(X_data_all, Xdash_data_all, library.F, options);
        end        
        
        
    case {'least-squares','least_squares','ls','nonlinear','standard'}
        
        %%% GATHER TOGETHER DATA
        
        % Loop over each dataset and count observations (for initialisation)
        N_obs = zeros(1,N_data);
        for k = 1:N_data
            N_obs(k) = length(data{k}.t);
        end

        % Now initialise data storage
        t_data_all = zeros(1, sum(N_obs));
        X_data_all = zeros(N_feat, sum(N_obs));
        
        % Loop over each replicate, adding to the total set of observations
        start_loc = 0;
        for k = 1:N_data

            % Define location in observation list for these replicates
            loc = start_loc + (1:N_obs(k));
            % Store the data in the master datalist in this location
            t_data_all(loc) = data{k}.t;
            X_data_all(:,loc) = data{k}.X;
            % Advance start location in preparation for next loop
            start_loc = start_loc + N_obs(k);
    
        end
        
        
        %%% USE THE NONLINEAR REPLICATOR FITTER TO SET COEFFICIENTS
        
        % This takes place in a separate function
        library.coeffs = nonlinearReplicatorFit(t_data_all, X_data_all, library, selection_type, options);
        
end