function dX = libraryRHS(X,F,K)
%
%   dX = libraryRHS(X,F,K)
%
% This function acts as an ODE right hand side function for use in solvers
% like ode15s, evaluating the ODE system as defined by an input set of RHS
% functions F, provided as a cell array. The vector K specifies the
% coefficients of each function.

% Evaluate the first function to set the size of dX
dX = K(1) * F{1}(X);

% Add on the contributions from the remaining library functions
for k = 2:length(K)
    dX = dX + K(k) * F{k}(X);
end