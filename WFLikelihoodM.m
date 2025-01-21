function loglike = WFLikelihoodM(data, problem)

% Read out the problem details
N = problem.N_pop;
F = problem.fitness;

% Ensure time starts at zero
data.t = data.t - data.t(1);

% Convert all frequency data into precise counts of individuals
data.X = round( N * data.X );

% Number of generations to run will be those contained within the data
N_gen = max(data.t);

% Store the number of states
Nstates = 1;
Nfeat = size(data.X,1);
for m = 1:Nfeat-1
    Nstates = Nstates * (N+m)/m;
end

% Loop over all the different states and store their frequency values
looping = true;
statelist = zeros(Nstates,N);
state = ones(1,N);
statelist(1,:) = state;
loc = 1;
uploc = N;
while looping
        
    % If we cannot increment at the current location, move backwards/exit
    if state(uploc) == Nfeat
        if uploc > 1
            uploc = uploc - 1;
        else
            looping = false;
        end
        
    % If we can increment, store
    else
        
        % Update storage location
        loc = loc + 1;
        % Update state
        state(uploc) = state(uploc) + 1;
        state(uploc+1:end) = state(uploc);
        % Add current state to list of states
        statelist(loc,:) = state;
        % Reset update location to end of vector
        uploc = N; 
    end
    
end

% Store all the different counts and frequency values for each state
Cstates = zeros(Nfeat,Nstates);
for m = 1:Nfeat
    Cstates(m,:) = sum( statelist == m, 2 );
end
Xstates = Cstates / N;
            
% Calculate the fitness values for all possible states
for k = 1:Nstates
    Fstates(:,k) = F( Xstates(:,k) );
end

% Selection-driven proportions are given by normalise(x .* f), so calculate
% all of these values beforehand
if problem.selection_type == 1
    Pstates = Xstates .* Fstates + (1 - sum(Xstates .* Fstates,1)) .* Xstates;
else
    Pstates = Xstates .* Fstates ./ ( sum(Xstates .* Fstates,1) );
end

% Calculate the log factorials
logFactK = zeros(1,N+1);
for k = 1:N+1
    
    % Use Stirling's formula if the value is too big
    if k-1 < 125
        logFactK(k) = log(factorial(k-1));
    else
        logFactK(k) = (k-1) * log(k-1) - k + 1 + 0.5*log(2*pi*(k-1));
    end
    
end

% Use these factorial values to pre-calculate log permutation counts
logperms = logFactK(end);
for m = 1:Nfeat
    logperms = logperms - logFactK(Cstates(m,:)+1);
end

%%% BUILD THE TRANSITION MATRIX OF LOG PROBABILITIES

% %%% VECTORISED CONSTRUCTION IF MATRIX NOT TOO BIG
if Nstates < 10000
 
     % Each column begins as log number of permutations for that state
     logQ = repmat(logperms',1,Nstates);
     % Then, loop over each feature and add its probabilistic contribution
     for m = 1:Nfeat
         % Find non-zero contributions
         nz = Cstates(m,:) > 0;
         % Add the contribution of the term  log( p_m ^ i )
         logQ(nz,:) = logQ(nz,:) + Cstates(m,nz)' .* log( Pstates(m,:) );
     end
         
     % Convert log probabilities to actual probabilities in a sparse matrix
     Q = zeros(Nstates,Nstates);
     sig_ind = logQ >= -25;
     Q(sig_ind) = exp(logQ(sig_ind));
     
%%% ELEMENT-WISE CONSTRUCTION TO AVOID MEMORY ISSUES
else
     
    % Initialise sparse storage
    ivals = [];
    jvals = [];
    Qvals = [];
    
    % Do a manual loop where each element is included if sufficiently large
    for i = 1:Nstates
        
        % Read out current state, and list of non-zero species
        Chere = Cstates(:,i);
        Cnz = Chere>0;
        Chere = Chere(Cnz);
        
        % Find a list of all other states sufficiently close to this state
        dists = pdist2(Cstates(:,i)',Cstates');
        looplist = find(dists < 80);        
        
        % Cut down the list of probabilities to only those relevant
        Phere = Pstates(Cnz,looplist);
        
        % Prepare vectorised calculation of log probabilities
        logE = logperms(i) + sum(Chere .* log( Phere ), 1);
        
        % Cut down to only significant probabilities
        sig = (logE >= -20);
        num_sig = sum(sig);
        
        % Add these significant values to masterlist
        ivals(end+1:end+num_sig) = i;
        jvals(end+1:end+num_sig) = looplist(sig);
        Qvals(end+1:end+num_sig) = exp(logE(sig));
        
    end
    
    % Now store in the sparse matrix
    Q = sparse(ivals,jvals,Qvals, Nstates, Nstates);
        
 end
            
                
            

% Loop through the data and find what state each datapoint corresponds to
for k = 1:size(data.X,2)
    [~,data.S(k)] = ismember(data.X(:,k)',Cstates','rows');
end

% Extract first datapoint as starting point
Pmat = zeros(Nstates,1);
Pmat(data.S(1)) = 1;

% Repeatedly iterate the transition matrix to calculate updated
% probabilities
for m = 1:N_gen
    Pmat(:,m+1) = Q * Pmat(:,m);
end

% Cut down the probability matrix to only those generations for which a 
% result was observed in the data
Pmat = Pmat(:,data.t + 1);

% Log-likelihood is the sum of log probabilities of being in each precise
% state observed in the data
loglike = 0;
for k = 1:length(data.t)
    loglike = loglike + log( Pmat(data.S(k),data.t(k)+1) );
end