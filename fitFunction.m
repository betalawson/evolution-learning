function [F, Fdash] = fitFunction(x_data, y_data, options)
%
%    [F, Fdash] = fitFunction(x, y)
%    [F, Fdash] = fitFunction(x, y, options)
%
% This function takes as input a set of points (x,y), each specified by
% vectors, and attempts to fit a function to that data using the requested
% method. If no method is specified, Gaussian process regression is used.
%
% -Inputs-
%
%              x: An Ndata by 1 vector of x observations
%
%              y: An Ndata by 1 vector of y observations
%             
%        options: A struct with parameters that the user wishes not to be
%                 default, such as which method to use for estimating the
%                 derivative, parameters of that method, and whether or not
%                 to visualise the fitting results
%
%
% -Outputs-
%
%         F: A function of form F(X) that outputs a predicted y value for 
%            an input x value
%     Fdash: A function of form F(X) that outputs an estimate of the
%            derivative, calculated using F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% USER SPECIFIED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
options.N_pts = 1000;                % Number of points used in function plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert all data to columns
x_data = x_data(:);
y_data = y_data(:);

% Only visualise if specifically requested to
if nargin < 3
    options = imbueELdefaults([],1);
end

% Fit the function according to the requested method
switch lower(options.deriv_method)
    
    case {'gp','gpr'}
        
        % Fit a Gaussian process to the data using MATLAB's built-in
        GP = fitrgp(x_data,y_data,'BasisFunction','Constant');
        % Create a function that evaluates the fitted GP for F(x)
        F = @(x) reshape(predict(GP,x(:)),size(x));
        % Create a function that evaluates F'(x) via finite differencing
        Fdash = @(x) (F(x + sqrt(eps)) - F(x - sqrt(eps))) / (2*sqrt(eps));

    case {'gp-careful'}
        
        % Fit a Gaussian process to the data using MATLAB's built-in
        GP = fitrgp(x_data,y_data,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('Verbose',0,'ShowPlots',false));
        % Create a function that evaluates the fitted GP for F(x)
        F = @(x) reshape(predict(GP,x(:)),size(x));
        % Create a function that evaluates F'(x) via finite differencing
        Fdash = @(x) (F(x + sqrt(eps)) - F(x - sqrt(eps))) / (2*sqrt(eps));
   
    case {'local-polynomial','poly','local'}
                
        % Use the subfunction defined below to evaluate the derivative
        F = @(x) localPolynomialFit(x, x_data, y_data, options.kernel_bandwidth, options.polynomial_order, options.polyfit_regularisation, 0);
        Fdash = @(x) localPolynomialFit(x, x_data, y_data, options.kernel_bandwidth, options.polynomial_order, options.polyfit_regularisation, 1);
        
end

% Plot the predictions if requested
if options.visualise_fit
       
    % Prepare figure
    figure('units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    % Calculate and plot functional fit
    xpts = linspace(min(x_data), max(x_data),options.N_pts)';
    Ffit = F(xpts);
    plot(xpts, Ffit, 'r', 'LineWidth', 2.5);
    % Plot data over the top
    plot(x_data, y_data, 'k.', 'MarkerSize', 30);
    
end




function out = localPolynomialFit(x, x_data, y_data, h, p, lambda, est_deriv)
% This subfunction solves the least squares problem associated with "local
% polynomial regression"

% If provided with a null bandwidth, attempt to automatically set it using
% point separation
if h == -1
    h = mean(diff(x_data))*10;
end

% Repeat the below process for each input point (not vectorised ouch)
[N1,N2] = size(x);
out = zeros(N1,N2);
for n1 = 1:N1
    for n2 = 1:N2
        
        % Grab out the current evaluation point
        x0 = x(n1,n2);
        % Generate weights for each datapoint based on tricube kernel
        w = 70/81 * (1 - (abs((x_data - x0)/h)).^3 ).^3 / h;
        % Grab out only those points with w > 0
        valid = (w > 0);
        N_valid = sum(valid);
        % Build weight matrix
        W = diag( w(valid) );
        
        % Prepare storage for design matrix
        X = zeros(N_valid, p+1);
        
        % Build the reduced regression problem for this weight template
        X(:,1) = 1;
        for k = 1:p
            X(:,k+1) = (x_data(valid) - x0).^k;
        end
        
        % Solve the least squares problem to get the estimator
        beta = ( X' * W * X + lambda * eye(p+1) ) \ ( X' * W * y_data(valid) );
        % The function estimate is the sum up to requested order of terms
        out(n1,n2) = factorial(est_deriv) * beta(est_deriv+1);
        
    end
end