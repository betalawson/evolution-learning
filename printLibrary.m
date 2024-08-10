function printLibrary(K,texts,threshold)
% 
%     printLibrary(K, texts)
%     printLibrary(K, texts, threshold)
%
% This function takes texts corresponding to a set of library functions, as
% will have been constructed by calling a library constructor, and the
% coefficients for those library functions, and outputs these to the user
% to inform them of the system found by SINDy fitting.
%
% If a value for 'threshold' is provided, only terms for which the
% coefficient exceeds that threshold will be printed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default to a tiny threshold if none specified
if nargin < 3
    threshold = 10*eps;
end

% Create space to display the result
fprintf('\n\n- FITTED FUNCTION LIBRARY IS: \n|\n');

% Loop over all library functions with a non-zero coefficient, and output
% the text corresponding to that function
for i = find( abs(K) > threshold )'
    
    % Print out the reaction and its rate
    fprintf('| Value %g for term: %s\n',K(i),texts{i});
    
end

% Ending space
fprintf('|\n-----------------------------\n\n');

end