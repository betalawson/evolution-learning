function dX = replicatorRHS(X,F,selection_type)
%
%     dX = replicatorRHS(X,F,U)
%     dX = replicatorRHS(X,F,U,selection_type)
%
% This function acts as an ODE right hand side function for use in solvers
% like ode15s, evaluating the replicator equation with fitness function F,
% mutation matrix U. The user can optionally specify the type of selection
% to use, which changes the implicit normalisation and affects the rate of
% the dynamics.

% Calculate the fitness of each subspecies (allele value or phenotype)
f = F(X);

% Use these fitness values to define the current replicator update
if (nargin < 3) || (selection_type == 1)
    dX = X .* ( f - (f' * X) );
else
    dX = X .* ( f / (X' * f) - 1 );
end

% Use final row of dX to instead lock the solution to probability simplex
dX(end) = 1 - sum(X);    % Associated with a zero mass (algebraic)