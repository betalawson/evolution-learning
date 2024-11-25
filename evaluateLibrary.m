function val = evaluateLibrary(X,varargin)
%
%   val = evaluateLibrary(X,library)
%   val = evaluateLibrary(X,F,K)
%
% This function simply evaluates the combination of function pieces in the
% library, with coefficients as specified, at the point X. If provided a
% library object, the fields library.F and library.K define the functions,
% and their coefficients, respectively. Alternatively, F and K can be
% provided as separate inputs.
%
% This function can be used to construct ODE right hand side functions, or
% simply to evaluate the library. If one is working with a polynomial 
% library of order 2 or less, conversion to matrix form is likely to be 
% more efficient.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retrieve functions and coefficients from the input data
if nargin == 2
    F = varargin{1}.F; 
    K = varargin{1}.coeffs;
else
    F = varargin{1};
    K = varargin{2};
end
    
% Evaluate the first function to set the size of dX
val = K(1) * F{1}(X);

% Add on the contributions from the remaining library functions
for k = 2:length(K)
    val = val + K(k) * F{k}(X);
end