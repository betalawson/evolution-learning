function [A, Q] = libraryToPayoff(library,N_feat)
%
%    A = libraryToPayoff(library,N_feat)
%
% This function takes an input polynomial library, with coefficients 
% attached, and converts the library into the equivalent payoff matrix. If
% the second output Q is requested, the function also interrogates the
% library for quadratic terms, and builds the matrix associated with those.
% These matrices can be used to convert the modular function libraries back
% into something that can be stored and simulated more readily.
%
% NOTE: This function relies upon default naming for the dependent
% variables (x1, x2, ...). Specifying custom names or altering names of 
% library functions will break this.
%
%  INPUTS:
%
%        library - A struct with field 'texts', that contains a cell array
%                  holding descriptions of each individual library function
%                  and a field 'K' that gives the learned coefficient for
%                  each function
%
%         N_feat - The number of features in this problem. Providing this
%                  simply avoids needing to interrogate the library to try
%                  to learn how many features are present (which is not yet
%                  implemented!)
%
%
%  OUTPUTS:
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out the library components
texts = library.texts;
K = library.coeffs;

% Prepare storage
A = zeros(N_feat);                           % Linear (payoff), NxN matrix
Q = zeros(N_feat, N_feat*(N_feat+1)/2);      % Quadratic, tri(N)xN matrix      

% Loop over all strings in the library, checking what type of term they are
for k = 1:length(texts)
    
    % Check if this is a constant term
    present = strfind(texts{k},'Constant');
    if ~isempty(present)
        
        % Check which variable it applies to
        xloc = strfind(texts{k},'x');
        colonloc = strfind(texts{k},':');
        var = str2double( texts{k}(xloc(1)+1:colonloc-1) );
        
        % A constant contributes to all payoff matrix elements in that row
        A(var,:) = A(var,:) + K(k);
        
    end
    
    % Check if this is a linear term
    present = strfind(texts{k},'Linear');
    if ~isempty(present)
        
        % Check which variable it applies to 
        xloc = strfind(texts{k},'x');
        colonloc = strfind(texts{k},':');
        var1 = str2double( texts{k}(xloc(1)+1:colonloc-1) );
        
        % Check which variable it depends on
        var2 = str2double( texts{k}(xloc(2)+1:end) );
        
        % Store in corresponding payoff matrix location
        A(var1,var2) = A(var1,var2) + K(k);
        
    end
    
    % Check if this is a quadratic term
    present = strfind(texts{k},'Quadratic');
    if ~isempty(present)
       
        % Check which variable it applies to 
        xloc = strfind(texts{k},'x');
        colonloc = strfind(texts{k},':');
        var1 = str2double( texts{k}(xloc(1)+1:colonloc-1) );
        
        % Check which variables it depends on
        byloc = strfind(texts{k},'by');
        var2a = str2double( texts{k}(xloc(2)+1:byloc-2) );
        var2b = str2double( texts{k}(xloc(3)+1:end) );
        
        % Convert this variable pair into its column in the quadratic
        % matrix (this is a hand-derived formula)
        col = var2b - N_feat + var2a * (N_feat - 0.5*(var2a-1));
        
        % Add to the quadratic terms matrix, in row given by which feature
        % this fitness applies to
        Q(var1,col) = Q(var1,col) + K(k);
        
    end
    
    % Check if this is a symmetric payoff term
    present = strfind(texts{k},'Symmetric');
    if ~isempty(present)
       
        % First variable - number after x and before 'and'
        xloc = strfind(texts{k},'x');
        andloc = strfind(texts{k},' and');
        var1 = str2double( texts{k}(xloc(2)+1:andloc-1) );
        
        % Check which variable it depends on
        var2 = str2double( texts{k}(xloc(3)+1:end) );
        
        % Store in corresponding payoff matrix location(s)
        A(var1,var2) = A(var1,var2) + K(k);
        if var1 ~= var2
            A(var2,var1) = A(var2,var1) + K(k);
        end
        
    end
    
end