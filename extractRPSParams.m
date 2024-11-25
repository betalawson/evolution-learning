function params = extractRPSParams(A,type)
%
% This function extracts the parameters for a rock paper scissors problem
% from the provided payoff matrix "A", depending on the type of selection
% (which affects which transforms can be applied invariantly to the
% matrix). For Type I selection, the diagonal can always be set to a value
% of unity, which is the chosen value for all baseline feature fitnesses.
% For type II, the average scaling of all diagonal elements is applied 
% using the geometric average.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check dimensions of matrix
N = size(A,1);

% Type I selection: set all diagonals to unity by add/subtract a constant
% applied to each column in turn
if type == 1
    for c = 1:N
        A(:,c) = A(:,c) - A(c,c) + 1;
    end
    
% Type II selection: take the geometric average of the diagonals
else
   
    % Take the average scaling of the diagonal elements
    Aii_mean = geomean(diag(A));
    % Scale the matrix by this to make the diagonals approximately unity
    A = A / Aii_mean;
    
end

% Subtract from each row the diagonal element in order to get zeros on the
% diagonals
f0 = zeros(1,3);
for r = 1:N
    
    % Store the subtracted constant as baseline fitness
    f0(r) = A(r,r);
    % Actually subtract the constant
    A(r,:) = A(r,:) - A(r,r);
    
end

% Now, store in a struct the fitness base parameters
for i = 1:N
    params.(['f_',num2str(i)]) = f0(i);
end

% Store the remaining off-diagonal matrix elements
for i = 1:N
    for j = 1:N     
        % Off-diagonals only
        if i~=j
            params.(['a_',num2str(i),num2str(j)]) = A(i,j);
        end
    end
end