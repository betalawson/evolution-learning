function [F, Fdash] = fitSimplexConstrainedFunction(data, options)
% This function creates a fit, F(t), to the input data, which should be a
% struct containing a series of time points in field 't', and a series of
% observations in field 'X' that correspond to those timepoints. The field
% 'X' should be provided as a N_feat x N_obs matrix

%%%%%%%%%%%%%%%%%%%%%%%% USER SPECIFIED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
options.N_pts = 1000;                % Number of points used in function plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in default options if they weren't provided
if nargin < 2
    options = addDefaultOptions([],1);
else
    options = addDefaultOptions(options,1);
end

% Use the first point in the provided data to read out number of features
N_feat = size(data.X,1);

% Loop over each feature and fit a function to its data
f = cell(1,N_feat);
for n = 1:N_feat
    
    % Generate fit and return it as a function handle f(t)
    f{n} = fitFunction(data.t, data.X(n,:), options);

end

% Create a function that evaluates all of the individual functions, and
% projects the result onto the probability simplex
F = @(t) projectOntoSimplex( evaluateFunctions(t,f) );

% Create a corresponding function that approximates the derivative of F
h = 1e-5;
Fdash = @(t) ( F(t + h) - F(t - h) ) / (2 * h);

% If requested, plot the fit (used to check fit performance)
if options.visualise_fit

    % Generate the fit at lots of points across the data range
    tvec = linspace(0,max(data.t),options.N_pts);
    fitX = F(tvec);
    
    % Loop over each feature separately for the plot
    for n = 1:N_feat
    
        % Prepare figure
        figure('units','normalized','OuterPosition',[0 0 1 1]); hold on;
        
        % Plot the data and the fit for this feature
        plot(tvec, fitX(n,:), 'r', 'LineWidth', 2);
        plot(data.t, data.X(n,:), 'k.', 'MarkerSize', 30);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = evaluateFunctions(inputs,funs)
% This subfunction simply evaluates a cell array of scalar-valued functions
% at the input points and stores the result in one big matrix (this is
% essentially a stripped down version of cellfun)

vals = zeros(length(funs), length(inputs));
for n = 1:length(funs) 
    vals(n,:) = funs{n}(inputs);
end