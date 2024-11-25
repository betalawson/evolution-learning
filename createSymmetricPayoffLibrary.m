function library = createSymmetricPayoffLibrary(N_feat, names, interactions)
%
%    library = createSymmetricPayoffLibrary(N)
%    library = createSymmetricPayoffLibrary(N, names)
%    library = createSymmetricPayoffLibrary(N, [], interactions) 
%    library = createSymmetricPayoffLibrary(N, names, interactions)
%
% This function creates a set of linear frequency-dependence terms that
% correspond directly to a payoff matrix, with a condition of symmetry
% enforced. This is useful for payoff matrices corresponding to diploid
% organisms where individual alleles have associated fitness advantage and
% possible dominance in heterozygotes.
%
% The user can optionally supply the names of the features, and a matrix
% that specifies which features are allowed to interact in the payoff
% matrix
%
% -Inputs-
%
%              N: Number of species in the system
%
%              m: The maximum order of polynomials to use in the library
%                 *** Currently only values m = {0,1,2,3} are supported ***
%                 
%          names: A cell array of text naming the N features (optional)
%
%   interactions: An NxN logical matrix specifying which species are
%                 allowed to interact in mixed polynomial terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INITIAL PREPARATION

% If no names provided, assign generic names
if nargin < 2 || isempty(names)
    for k = 1:N_feat
        names{k} = ['x',num2str(k)];
    end
end

% If not provided, assume all interactions are allowed
if nargin < 3
    interactions = true(N_feat);
end

% Functions and their associated texts are stored in cell arrays. 
funs = {};
texts = {};

% Each species gets one fitness function, so loop over species
for k = 1:N_feat
    
    %%% All payoff matrix elements correspond to linear fitness pieces
    
    % Only include unique pairings, i.e. don't double count (1,2) and (2,1)
    for j = k:N_feat
        % Only include if species j can contribute to fitness of k
        if interactions(j,k)
            
            % Library functions are made up of two pieces on offdiagonals, corresponding to both elements that need to be modified equivalently for symmetry
            funs{end+1} = @(X) [zeros(k-1,size(X,2)); X(j,:); zeros(N_feat-k,size(X,2)) ] + (j~=k) * [zeros(j-1,size(X,2)); X(k,:); zeros(N_feat-j,size(X,2)) ];
            texts{end+1} = ['Symmetric payoff matrix interaction for ',names{k},' and ',names{j}];
        
        end
    end

end

% Store these values in the library object
library.F = funs;
library.texts = texts;

end