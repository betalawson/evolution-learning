function loglike = WFLikelihood2(data, problem)

% Read out the problem details
N = problem.N_pop;
F = problem.fitness;

% Ensure time starts at zero
data.t = data.t - data.t(1);

% Convert all frequency data into precise counts of individuals
data.X = round( N * data.X );

% Number of generations to run will be those contained within the data
N_gen = max(data.t);

% Check if this is m = 2 data as only that is currently implemented
if size(data.X,1) ~= 2
    error('Provided data does not have two features!');
end

% Store all the different frequency values for each state
Xstates = linspace(0,1,N+1);
Xstates(2,:) = 1 - Xstates(1,:);

% Calculate the fitness values for all possible states
for k = 1:N+1
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
        logFactK(k) = (k-1) * log(k-1) - k + 1 + 0.5*log(2*pi*(k-1)) + 1/(12*(k-1)) - 1 / (360*(k-1)^3);
    end
    
end

% Use these factorial values to evaluate log permutation counts
logperms = logFactK(end) - logFactK(1:N+1) - logFactK(N+1:-1:1);

% Build the full matrix of transition probabilities
logQ = logperms' + (0:N)' .* log( Pstates(1,2:end-1) ) + (N:-1:0)' .* log( Pstates(2,2:end-1) );
%Q = zeros(N+1,N-1);
%sig_ind = logQ >= -25;
%Q(sig_ind) = exp(logQ(sig_ind));
Q = exp(logQ);
% Append the first and last columns to the matrix (no transitions out of
% absorbing states)
Q = [ [1;zeros(N,1)], Q, [zeros(N,1);1] ];

% Extract first datapoint as starting point
Pmat = zeros(N+1,1);
Pmat(data.X(1,1)+1) = 1;

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
    loglike = loglike + log( Pmat(data.X(1,k)+1,k) );
end