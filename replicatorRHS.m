function dX = replicatorRHS(X,F,selection_type,DAE)
%
%     dX = replicatorRHS(X,F)
%     dX = replicatorRHS(X,F,selection_type)
%     dX = replicatorRHS(X,F,selection_type,DAE)
%
% This function acts as an ODE right hand side function for use in solvers
% like ode15s, evaluating the replicator equation with fitness function F,
% mutation matrix U. The user can optionally specify the type of selection
% to use, which changes the implicit normalisation and affects the rate of
% the dynamics.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the fitness of each subspecies (allele value or phenotype)
f = F(X);

% Use these fitness values to define the current replicator update
if (nargin < 3) || (selection_type == 1)
    dX = X .* ( f - (X' * f) );
else
    dX = X .* ( f / (X' * f) - 1 );
end

% If in DAE mode, use the final row of dX to instead lock the solution to
% the probability simplex (this row is associated with a zero mass). This
% is not default behaviour, must be specified
if nargin < 4
    DAE = false;
end
if DAE  
    dX(end) = 1 - sum(X);
end