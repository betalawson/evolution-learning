function F = fitnessRPS(N_feat, f_base, b, c)
%
%    F = fitnessRPS(N_feat, f_base, b, c)
%
% This function constructs a fitness function of form f = Ax, where the
% payoff matrix A corresponds to general rock paper scissors dynamics.
% The strength of each strategy against the strategy it beats is given by
% input vector "b", whereas the extent to which each strategy loses to the
% strategy that beats it is given by "c". Both are specified using positive
% values. If the user provides only a single value for "b" and/or "c", that
% value is used for all elements of the corresponding vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set b and c to vectors of appropriate size if scalar input provided
if isscalar(b)
    b = b * ones(N_feat,1);
end
if isscalar(c)
    c = c * ones(N_feat,1);
end

% Initialise the payoff matrix to the base fitness value
A = f_base * ones(N_feat);

% Create a cyclic function for matrix element references
ind1mod = @(x) mod(x-1,N_feat)+1;

% Loop through matrix elements, creating species that act negatively and
% positively on others
count = 0;
for i = 1:N_feat
    
    pos_spec = ind1mod(i+1);
    neg_spec = ind1mod(i+2);
    count = count + 1;
    
    A(i,pos_spec) = A(i,pos_spec) + b(count);
    A(i,neg_spec) = A(i,neg_spec) - c(count);

end

% Create the fitness function using the constructed payoff matrix
F = @(X) A * X;