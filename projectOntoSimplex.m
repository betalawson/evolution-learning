function X = projectOntoSimplex(X)
% This function projects an input column vector X onto the probabilty
% simplex, that is, it finds the closest (in terms of Euclidean distance) 
% vector X that satisfies sum(X) = 1 and all elements of X in [0,1]. 
% This algorithm is from Wang 2013 (citation given in associated paper)
%
% The function may also take an input matrix X, where each column is a Dx1
% vector to be projected.

% Read out the dimension of the input
D = size(X,1);

% Sort the elements of (each column of) X
u = sort(X,1,'descend');

% Calculate the running sum of elements used to decide where to truncate
rho_check = diag(1./(1:D)) * (1 - cumsum(u));

% Now just process each column independently for easy code
for c = 1:size(X,2)
    
    % Find the location just before this shift applied to current element
    % switches to negative
    rho = find( u(:,c) + rho_check(:,c) > 0, 1, 'last');

    % Get the shift corresponding with this location
    lambda = rho_check(rho,c);

    % Apply this shift to all elements of this column, but ensure positive
    X(:,c) = max( [X(:,c) + lambda, zeros(D,1)], [], 2 );
    
end