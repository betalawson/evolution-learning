function [funs, texts] = createPolynomialLibrary(N_feat, orders, names, interactions)
%
%    [funs, texts] = createPolynomialLibrary(N, m)
%    [funs, texts] = createPolynomialLibrary(N, m, names)
%    [funs, texts] = createPolynomialLibrary(N, m, [], interactions) 
%    [funs, texts] = createPolynomialLibrary(N, m, names, interactions)
%
% This function creates a set of polynomial library functions, composed of
% all polynomials up to a specified order that obey the specified structure
% of terms that are allowed to interact (or all polynomials if no
% interaction matrix is provided)
%
% -Inputs-
%
%              N: Number of features in the system
%
%              m: The maximum order of polynomials to use in the library
%                 *** Currently only values m = {0,1,2,3} are supported ***
%                 
%          names: A cell array of text naming the N species (optional)
%
%   interactions: An NxN logical matrix specifying which species are
%                 allowed to interact in mixed polynomial terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INITIAL PREPARATION

% If no names provided, assign generic names
if nargin < 3 || isempty(names)
    for k = 1:N_feat
        names{k} = ['x',num2str(k)];
    end
end

% If not provided, assume all interactions are allowed
if nargin < 4
    interactions = true(N_feat);
end

% Functions and their associated texts are stored in cell arrays. 
funs = {};
texts = {};

if max(orders) > 3
    warning('Polynomial orders greater than three are currently not supported! Will only learn polynomial orders up to 3.');
end

% Each species gets one fitness function, so loop over species
for k = 1:N_feat
    
    %%% Zeroth-order terms
    if ismember(0,orders)
        funs{end+1} = @(X) [zeros(k-1,size(X,2)); ones(1,size(X,2)); zeros(N_feat-k,size(X,2))];
        texts{end+1} = ['Fitness for ',names{k},': Constant'];
    end
    
    %%% First-order terms
    if ismember(1,orders)
        for j = 1:N_feat
            % Only include if species j can contribute to fitness of k
            if interactions(j,k)
                funs{end+1} = @(X) [zeros(k-1,size(X,2)); X(j,:); zeros(N_feat-k,size(X,2)) ];
                texts{end+1} = ['Fitness for ',names{k},': Linear term - ',names{j}];
            end
        end
    end
    
    %%% Second-order terms
    if ismember(2,orders)
        for i = 1:N_feat
            for j = i:N_feat
                % Only include if species i,j can contribute to fitness of k
                if interactions(i,k) && interactions(j,k) && interactions(i,j)
                    funs{end+1} = @(X) [zeros(k-1,size(X,2)); X(i,:) .* X(j,:); zeros(N_feat-k,size(X,2)) ];
                    texts{end+1} = ['Fitness for ',names{k},': Quadratic term - ',names{i},' by ',names{j}];
                end
            end
        end
    end
    
    %%% Third-order terms
    if ismember(3,orders)
        for i = 1:N_feat
            for j = i:N_feat
                for l = j:N_feat
                    % Only include if species i,j,l can contribute to fitness of k
                    if interactions(i,k) && interactions(j,k) && interactions(l,k) && interactions(i,j) && interactions(i,l) && interactions(j,l)
                        funs{end+1} = @(X) [zeros(k-1,size(X,2)); X(i,:) .* X(j,:) .* X(l,:); zeros(N_feat-k,size(X,2)) ];
                        texts{end+1} = ['Fitness for ',names{k},': Cubic term - ',names{i},' by ',names{j},' by ',names{l}];
                    end
                end
            end
        end
    end
    
end

end
