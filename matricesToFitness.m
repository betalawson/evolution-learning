function F = matricesToFitness(A,Q)
%
%    F = matricesToFitness(A)
%    F = matricesToFitness(A,Q)
%
% This function simply returns a single function, F, that represents the
% product of the payoff matrix A with the feature frequencies, and the
% quadratic matrix Q (if provided and non-zero) with the unique
% combinations of feature frequencies (second-order terms)
%
%  INPUTS:
%
%              A - A matrix of dimension N_feat x N_feat that represents
%                  the linear dependencies of the fitness of each feature,
%                  on the frequencies of each of the features. Note that
%                  constant fitnesses of features are also contained within
%                  each row of this matrix
%
%              Q - A matrix of dimension N_feat*(N_feat-1)/2 x N_feat that
%                  contains the quadratic dependencies of the fitness of
%                  each feature, on the frequencies of each unique product
%                  of features. Note that A and Q are not unique, as terms
%                  can be shifted from A to Q or vice versa in similar
%                  fashion to how constant terms can be brought into/out
%                  from A
%
%  OUTPUTS:
%
%              F - A function of form F(X) that evaluates the fitness of
%                  all features, given their current frequencies X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Use the input matrix to check the dimension
N = size(A,1);

% First, check if any quadratic dependencies exist
if nargin < 2
    quadratic = false;
elseif any( abs(Q) > 1e-12 )
    quadratic = true;
else
    quadratic = false;
end

% If the quadratic term is not present, skip those terms
if ~quadratic
    F = @(X) A * X;
    
% Otherwise, use the full expression
else
    
    %%% Pre-generate a list of unique pairs of features
    
    % Create matrices for values 1:N along each row, and along each column
    I = repmat(1:N,N,1);
    J = I';
    % Read out the full list of unique pairs from these
    first_feat = I(I<=J);
    second_feat = J(I<=J);
    
    % Now get products of the relevant matrices with linear and quadratic
    % combinations in X
    F = @(X) A * X + Q * (X(first_feat) .* X(second_feat));
   
end