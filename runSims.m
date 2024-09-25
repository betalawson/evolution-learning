function simResults = runSims( problem, ELoptions, param_fun, pop_sizes, N_rep )
%
% This function generates multiple (N_rep) copies of data from the
% Wright-Fisher model for the provided problem object 'problem', for each
% of the different population sizes provided in the input vector
% 'pop_sizes'. Then, for each of those generated datasets, both standard
% and gradient-matching equation learning are applied (with options in
% provided object ELoptions). Parameters are extracted from the learned
% dynamics using the provided function "param_fun", which should be a
% function that returns a structure with fields being the parameters, and
% takes as inputs (K,Ftexts,type):
%
%           K - library coefficients
%      Ftexts - library texts defining what each function is
%        type - selection type (value of 1 or 2)

% Specify a random seed for consistency - this should give the same
% trajectories for each different approach to match, also
rng(7);

% Extract number of population sizes to run over
N_pops = length(pop_sizes);


%%% REPEATED RUNS ACROSS WRIGHT-FISHER REALISATIONS

% Initialise storage - separate for both approaches (standard and
% gradient-matching) due to use of parfor
l2LS_mat = zeros(N_pops, 2, N_rep);
l2EL_mat = zeros(N_pops, 2, N_rep);
paramsLS_mat = cell(N_pops, 2, N_rep);
paramsEL_mat = cell(N_pops, 2, N_rep);
fLS_mat = cell(N_pops, 2, N_rep);
fEL_mat = cell(N_pops, 2, N_rep);

% Initialise total time structure
timing = struct('LS',0,'EL',0);

% Now, loop over the list of population sizes
for k = 1:N_pops
    
    % Loop over the two types of selection and generate results for each
    for type = 1:2
        
        % Prepare problem object using input problem but with current
        % selection type and population size
        Pobj = setfield(problem, 'selection_type', type);
        Pobj.N_pop = pop_sizes(k);
        
        % Initialise storage for timing
        LS_toc = zeros(1,N_rep);
        EL_toc = zeros(1,N_rep);
        
        % Repeatedly generate Wright-Fisher data, measure fit, get params
        parfor r = 1:N_rep
            
            % Generate the Wright-Fisher and replicator data
            traj = generateTrajectories(Pobj);
            
            % Convert to non-dimensionalised time for equation learning
            WF_nonD = struct('t',traj.WF.t / Pobj.t_gen, 'X', traj.WF.X);
            
            % Perform least-squares equation learning
            tic;
            [Flib, K, Ftexts] = LSevolutionLearning(WF_nonD, type, ELoptions);
            LS_toc(r) = toc;
            
            % Store fitness function and params for the learned dynamics
            F_here = @(X) libraryRHS(X, Flib, K);
            fLS_mat{k,type,r} = F_here;
            paramsLS_mat{k,type,r} = param_fun(K, Ftexts, type);
            
            % Compare learned fitness to the data
            learned_traj = generateTrajectories( setfield(Pobj,'fitness', @(X) libraryRHS(X, Flib, K)), 'rep' );
            if isequal(size(traj.rep.X(:)),size(learned_traj.rep.X(:))) 
                l2LS_mat(k,type,r) = rms( traj.rep.X(:) - learned_traj.rep.X(:) );
            % If integration still managed to fail despite adaptive
            % tolerance specification, store a NaN (and the corresponding
            % entry in f_mat can be checked manually later)
            else
                l2LS_mat(k,type,r) = NaN;
            end
            
            % Now use our gradient matching method
            tic;
            [~, Flib, K, Ftexts] = evolutionLearning(WF_nonD, type, ELoptions);
            EL_toc(r) = toc;
            
            % Store fitness function and params for the learned dynamics
            F_here = @(X) libraryRHS(X, Flib, K);
            fEL_mat{k,type,r} = F_here;
            paramsEL_mat{k,type,r} = param_fun(K, Ftexts, type);
            
            % Compare learned fitness to the data
            learned_traj = generateTrajectories( setfield(Pobj,'fitness', @(X) libraryRHS(X, Flib, K)), 'rep' );
            if isequal(size(traj.rep.X(:)),size(learned_traj.rep.X(:)))
                l2EL_mat(k,type,r) = rms( traj.rep.X(:) - learned_traj.rep.X(:) );    
            % If integration still managed to fail despite adaptive
            % tolerance specification, store a NaN (and the corresponding
            % entry in f_mat can be checked manually later)
            else
                l2EL_mat(k,type,r) = NaN;
            end
            
        end
        
        % Increment total time counters
        timing.LS = timing.LS + sum(LS_toc);
        timing.EL = timing.EL + sum(EL_toc);
        
    end
end

% Gather together the data for least-squares and gradient matching
l2_mat(:,1:2,:) = l2LS_mat;
l2_mat(:,3:4,:) = l2EL_mat;
f_mat(:,1:2,:) = fLS_mat;
f_mat(:,3:4,:) = fEL_mat;
params_mat(:,1:2,:) = paramsLS_mat;
params_mat(:,3:4,:) = paramsEL_mat;


%%% ONE RUN PER METHOD/TYPE FOR PERFECT DATA

% Initialise storage - separate for both approaches (standard and
% gradient-matching) due to use of parfor
l2_perfect = zeros(1, 4);
params_perfect = cell(1, 4);

% Loop over the two types of selection and generate results for each
for type = 1:2
        
    % Prepare problem object using input problem but with current
    % selection type and population size
    Pobj = setfield(problem, 'selection_type', type);
           
    % Generate the replicator data including perfectly correct derivatives
    traj = generateTrajectories(Pobj,'rep');
    X_true = traj.rep.X;
    Xdash_true = zeros(size(X_true));
    for j = 1:size(X_true,2)
        Xdash_true(:,j) = replicatorRHS(X_true(:,j),problem.fitness,type);
    end
            
    % Calculate the best-fit parameters 
    [Flib, K, Ftexts] = LSevolutionLearning(struct('t',traj.rep.t/problem.t_gen,'X',traj.rep.X), type, ELoptions);
    
    % Extract parameters
    params_perfect{type} = param_fun(K, Ftexts, type);
    
    % Find root mean square error for these dynamics
    learned_traj = generateTrajectories( setfield(Pobj,'fitness', @(X) libraryRHS(X, Flib, K)), 'rep' );
    l2_perfect(type) = rms( traj.rep.X(:) - learned_traj.rep.X(:) );
    
    % Now use our gradient matching method, with the derivative estimates already provided
    traj.rep.Xdash = Xdash_true;
    [~, Flib, K, Ftexts] = evolutionLearning(traj.rep, type, ELoptions);
    
    % Extract parameters
    params_perfect{type+2} = param_fun(K, Ftexts, type);
    
    % Find root mean square error for these dynamics
    learned_traj = generateTrajectories( setfield(Pobj,'fitness', @(X) libraryRHS(X, Flib, K)), 'rep' );
    l2_perfect(type+2) = rms( traj.rep.X(:) - learned_traj.rep.X(:) );
    
end

% Store all the results in a struct
simResults = struct('l2_mat',l2_mat,'l2_perfect',l2_perfect,'params_mat',{params_mat},'params_perfect',{params_perfect},'f_mat',{f_mat},'pop_sizes',pop_sizes,'timing',timing);

end