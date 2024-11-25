function F = convertToReplicatorLibrary(F)
% This function takes as input a series of library functions that make up
% the definition of fitness, and converts them into a series of library
% functions for the replicator equation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out how many functions are here
Nfuns = length(F);

% Now loop over all the functions, and replace each with its
% replicator-transformed equivalent
for k = 1:Nfuns
    F{k} = @(X) X .* ( F{k}(X) - sum( X .* F{k}(X) ) .* ones(size(X)) );
end

end