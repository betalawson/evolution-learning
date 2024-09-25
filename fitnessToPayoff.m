function payoff = fitnessToPayoff(K,texts,N_spec)
% This function combines together the coefficients for constants and linear
% terms in the set of library functions with coefficients in (K,texts) and
% converts them into the equivalent payoff matrix.
%
% NOTE: This function relies upon default naming. Specifying names for
% variables or altering names of library functions will break this.

% Prepare storage
payoff = zeros(N_spec);
const = zeros(N_spec,1);
C = cell(1,N_spec);
b = cell(1,N_spec);
for n = 1:N_spec
    C{n} = zeros( round(N_spec*(N_spec+1)/2), N_spec );
    b{n} = zeros( round(N_spec*(N_spec+1)/2), 1 );
end
q = zeros(N_spec,1);

% Loop over all strings in library
for k = 1:length(texts)
    
    % Check if this is a constant term
    present = strfind(texts{k},'Constant');
    if ~isempty(present)
        
        % Check which variable it applies to
        xloc = strfind(texts{k},'x');
        colonloc = strfind(texts{k},':');
        var = str2double( texts{k}(xloc(1)+1:colonloc-1) );
        
        % Store in constant vector storage
        const(var) = K(k);
        
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
        payoff(var1,var2) = K(k);
        
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
        
        % Add this row to the matrix system used to "undo" quadratic terms,
        % by including a 1 in all columns this term relates to and in the
        % RHS vector including the coefficient for this term
        q(var1) = q(var1) + 1;
        C{var1}(q(var1),var2a) = 1;
        C{var1}(q(var1),var2b) = 1;
        b{var1}(q(var1)) = K(k);
        
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
        payoff(var1,var2) = K(k);
        payoff(var2,var1) = K(k);
        
    end
    
end

% Combine constants with payoff matrix (adding baseline fitness for a
% feature is equivalent to adding that value to all elements in its row)
payoff = const + payoff;

% If any quadratic terms exist, include these as well by solving the linear
% systems constructed to find optimal shifts to shift some portion of the
% quadratic terms to their linear equivalents
%if any(q > 0)
%    for n = 1:N_spec
%        payoff(n,:) = payoff(n,:) + ( C{n} \ b{n} )';
%    end
%end