function f = fitFunction(x_data, y_data, options)
%
%    f = fitFunction(x, y)
%    f = fitFunction(x, y, options)
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
%         f: A function of form F(X) that outputs a predicted y value for 
%            an input x value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert all data to columns
x_data = x_data(:);
y_data = y_data(:);

% Read in default options if they weren't provided
if nargin < 3
    options = imbueELDefaults([],1);
else
    options = imbueELDefaults(options,1);
end

% Fit the function according to the requested method
switch lower(options.deriv_method)
    
    case {'gp','gpr'}
        
        % Fit a Gaussian process to the data using MATLAB's built-in
        GP = fitrgp(x_data,y_data);
        % Create a function that evaluates the fitted GP for F(x)
        f = @(x) reshape(predict(GP,x(:)),size(x));
        
    case {'gp-careful'}
        
        % Fit a Gaussian process to the data using MATLAB's built-in
        GP = fitrgp(x_data,y_data,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('Verbose',0,'ShowPlots',false));
        % Create a function that evaluates the fitted GP for F(x)
        f = @(x) reshape(predict(GP,x(:)),size(x));
   
    case {'local-polynomial','poly','local'}
                
        % Use the subfunction defined below to evaluate the derivative
        f = @(x) localPolynomialFit(x, x_data, y_data, options.kernel_bandwidth, options.polynomial_order, options.polyfit_regularisation, 0);
        
end


function out = localPolynomialFit(x, x_data, y_data, h, p, lambda, est_deriv)
% This subfunction solves the least squares problem associated with "local
% polynomial regression", and returns an estimate of the derivative of
% order equal to the input variable "est_deriv".

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
        %w = 70/81 * (1 - (abs((x_data - x0)/h)).^3 ).^3 / h;
        w = 1/h * exp( -(x_data - x0).^2 / (2*h^2) );
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