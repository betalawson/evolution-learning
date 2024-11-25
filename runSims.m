function simResults = runSims( problem, ELoptions, param_fun, pop_sizes, N_rep, sim_gen )
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
obsl2LS_mat = zeros(N_pops, 2, N_rep);
obsl2EL_mat = zeros(N_pops, 2, N_rep);
predl2LS_mat = zeros(N_pops, 2, N_rep);
predl2EL_mat = zeros(N_pops, 2, N_rep);
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
            WF_nonDs = struct('t',traj.WF.t / Pobj.t_gen, 'X', traj.WF.X);
            
            % Perform least-squares equation learning
            tic;
            library = evolutionLearning(WF_nonDs, type, setfield(ELoptions,'regression_type','ls'));
            LS_toc(r) = toc;
            
            % Store fitness function and params for the learned dynamics
            F_here = @(X) evaluateLibrary(X, library);
            fLS_mat{k,type,r} = F_here;
            paramsLS_mat{k,type,r} = param_fun(library, type);
            
            % Calculate L2 error between learned and true replicator
            obsl2LS_mat(k,type,r) = calculateTrajectoryError(Pobj, F_here, Pobj.X0);
            predl2LS_mat(k,type,r) = calculateTrajectoryError(setfield(Pobj,'N_gen',sim_gen), F_here, Pobj.X0);
            
            % Now use our gradient matching method
            tic;
            library = evolutionLearning(WF_nonDs, type, setfield(ELoptions,'regression_type','gm'));
            EL_toc(r) = toc;
            
            % Store fitness function and params for the learned dynamics
            F_here = @(X) evaluateLibrary(X, library);
            fEL_mat{k,type,r} = F_here;
            paramsEL_mat{k,type,r} = param_fun(library, type);
            
            % Calculate L2 error between learned and true replicator
            obsl2EL_mat(k,type,r) = calculateTrajectoryError(Pobj, F_here, Pobj.X0);
            predl2EL_mat(k,type,r) = calculateTrajectoryError(setfield(Pobj,'N_gen',sim_gen), F_here, Pobj.X0);
            
        end
        
        % Increment total time counters
        timing.LS = timing.LS + sum(LS_toc);
        timing.EL = timing.EL + sum(EL_toc);
        
    end
end

% Gather together the data for least-squares and gradient matching
obsl2_mat(:,1:2,:) = obsl2LS_mat;
obsl2_mat(:,3:4,:) = obsl2EL_mat;
predl2_mat(:,1:2,:) = predl2LS_mat;
predl2_mat(:,3:4,:) = predl2EL_mat;
f_mat(:,1:2,:) = fLS_mat;
f_mat(:,3:4,:) = fEL_mat;
params_mat(:,1:2,:) = paramsLS_mat;
params_mat(:,3:4,:) = paramsEL_mat;


%%% ONE RUN PER METHOD/TYPE FOR PERFECT DATA

% Initialise storage - separate for both approaches (standard and
% gradient-matching) due to use of parfor
obsl2_perfect = zeros(1, 4);
predl2_perfect = zeros(1, 4);
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
    library = evolutionLearning(struct('t',traj.rep.t/problem.t_gen,'X',traj.rep.X), type, setfield(ELoptions,'regression_type','ls'));
    
    % Extract parameters and fitness function
    params_perfect{type} = param_fun(library, type);
    F_here = @(X) evaluateLibrary(X, library);
    
    % Find root mean square error for these dynamics
    obsl2_perfect(type) = calculateTrajectoryError(Pobj, F_here, Pobj.X0);
    predl2_perfect(type) = calculateTrajectoryError(setfield(Pobj,'N_gen',sim_gen), F_here, Pobj.X0);
        
    % Now use our gradient matching method, with the derivative estimates already provided
    traj.rep.Xdash = Xdash_true;
    library = evolutionLearning(traj.rep, type, setfield(ELoptions,'regression_type','grad'));
    
    % Extract parameters and fitness function
    params_perfect{type+2} = param_fun(library, type);
    F_here = @(X) evaluateLibrary(X, library);
    
    % Find root mean square error for these dynamics
    obsl2_perfect(type+2) = calculateTrajectoryError(Pobj, F_here, Pobj.X0);
    predl2_perfect(type+2) = calculateTrajectoryError(setfield(Pobj,'N_gen',sim_gen), F_here, Pobj.X0);
    
end

% Store all the results in a struct
simResults = struct('obsl2_mat',obsl2_mat,'obsl2_perfect',obsl2_perfect,'predl2_mat',predl2_mat,'predl2_perfect',predl2_perfect,'params_mat',{params_mat},'params_perfect',{params_perfect},'f_mat',{f_mat},'pop_sizes',pop_sizes,'timing',timing);

end

function l2 = calculateTrajectoryError(Pobj, F, X0list)
% This function generates replicator trajectories for the provided problem
% object, and for the provided object but with the specified fitness
% function F, and evaluates the average l2 error between them across all of
% the initial conditions provided in X0list

% Read out the number of ICs given (number of column vectors in X0list)
N_ICs = size(X0list,2);

% Loop over all trajectories and evaluate each's contribution to mean error
l2 = 0;
for k = 1:N_ICs
               
    % Create a problem object for this initial condition
    Phere = setfield(Pobj, 'X0', X0list(:,k));
    
    % Create trajectory for true selective dynamics
    true_traj = generateTrajectories( Phere, 'rep' );
    % Create trajectory for learned dynamics
    learned_traj = generateTrajectories( setfield(Phere, 'fitness', F), 'rep' );
                
    % Calculate error - if integration failed despite adaptive tolerance
    % specification, store a NaN
    if isequal(size(true_traj.rep.X(:)),size(learned_traj.rep.X(:))) 
        l2 = l2 + rms( true_traj.rep.X(:) - learned_traj.rep.X(:) ) / N_ICs;
    else
        l2 = NaN;
    end
                
end

end