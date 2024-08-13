function F = fitnessAlleles(s, H, s0)
%
%    F = fitnessAlleles(s, h)
%    F = fitnessAlleles(s, H)
%
% This function constructs a fitness function of form f = Ax, where the
% payoff matrix A corresponds to alleles in diploid organisms. The input
% vector 's' specifies the level of selection for the different values for
% the allele, while the matrix H is a set of values [0,1] that specifies
% the level of dominance in each pairing (for example, the element H(1,2)
% specifies the extent to which allele a1 dominates over a2 in the a1a2
% individual. Necessarily, H(i,j) = 1 - H(j,i). Or if a scalar 'h' is
% provided, this level of dominance is used for every element.

% Assume a value of s0 = 1 if none given
if nargin < 3
    s0 = 1;
end

% Use the provided selection strength vector to get number of species
N_spec = length(s);

% Ensure that H is of valid structure, using its lower diagonal
if ~isequal(size(H),[1 1])
    H = tril(H,-1) + triu((1 - tril(H,-1)'),1);
    
% If only a single scalar provided, make a 
else
    H = H * ones(N_spec);
    H = H - diag(diag(H));
end

% Build the whole payoff matrix
A = zeros(N_spec);
for i = 1:N_spec
    
    % Diagonal elements: s0 + s_i
    A(i,i) = s0 + s(i);
    
    % Off-diagonal elements: s0 + h(i,j) s(i) + h(j,i) s(j)
    for j = i+1:N_spec
        A(i,j) = s0 + H(j,i) * s(i) + H(i,j) * s(j);
        A(j,i) = s0 + H(j,i) * s(i) + H(i,j) * s(j);
    end
    
end

% Create the fitness function using the constructed payoff matrix
F = @(X) A * X;