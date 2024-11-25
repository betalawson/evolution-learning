function params = extractBialleleParams(A,type)
%
% This function extracts the parameters (s,h) from the payoff matrix that
% was provided as input, using the invariant rules allowed by this type
% of selection are used to normalise and (if possible) symmeterise the 
% matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalise and symmeterise the matrix
if type == 1
    A = A - A(2,2) + 1;
    A(:,1) = A(:,1) - A(2,1) + A(1,2);
else
    A = A / A(2,2);
    A = (A+A')/2;             % Forcibly symmeterise for type II
end
            
% Extract the values of s and h 
params.s = A(1,1) - 1;
params.hs = A(1,2) - 1;
