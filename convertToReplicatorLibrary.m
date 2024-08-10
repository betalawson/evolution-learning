function replicator_funs = convertToReplicatorLibrary(funs)
% This function takes as input a series of library functions that make up
% the definition of fitness, and converts them into a series of library
% functions for the replicator equation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out how many functions are here
Nfuns = length(funs);

% Now loop over all the functions, and replace each with its
% replicator-transformed equivalent
replicator_funs = cell(1,Nfuns);
for k = 1:Nfuns
    replicator_funs{k} = @(X) X .* ( funs{k}(X) - sum( X .* funs{k}(X) ) .* ones(size(X)) );
end

end