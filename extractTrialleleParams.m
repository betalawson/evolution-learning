function params = extractTrialleleParams(A,type)
%
% This function extracts the parameters (s1, s2, h12, h13, h23) from the 
% payoff matrix that was provided as input, using the invariant rules
% allowed by this type of selection to normalise the matrix. Note that the
% matrix is treated as symmetric as an assumption.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalise the matrix using its bottom-right element
if type == 1
    A = A - A(3,3) + 1;
else
    A = A / A(3,3);
end
            
% Extract the values of s and h 
params.s1 = A(1,1) - 1;
params.s2 = A(2,2) - 1;
params.s12 = A(1,2) - 1 - params.s2;
params.s13 = A(1,3) - 1;
params.s23 = A(2,3) - 1;
