function trajectories = generateTrajectories( problems, gen_traj )
% This function generates both a Wright-Fisher trajectory and the
% replicator trajectory for the problem (or set of problems if given a cell
% array) that are provided as input. Part of these problem objects should
% also specify how many generations to run the dynamics for, and how
% frequently these dynamics (an integer number of generations) should be 
% observed.
%
% The user can optionally specify a list of types of trajectories to
% generate, from the list 'wf', 'rep', 'both'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 N_simpts = 10001;   % Number of points to use for replicator trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, check if a cell array was provided
if iscell(problems)
    cell_input = true;
    N_prob = length(problems);
% If an array wasn't provided, create a length-one array for consistency
else
    cell_input = false;
    N_prob = 1;
    problems = {problems};
end

% If a list of trajectories to generate was not provided, assume all
if nargin < 2
    gen_traj = 'both';
end

% Loop over each problem, creating the trajectory object for each
trajectories = cell(N_prob,1);
for k = 1:N_prob
   
    % Read out all the current problem's settings
    problem = problems{k};
    N_feat = problem.N_feat;
    fitness = problem.fitness;
    X0 = problem.X0;
    N_pop = problem.N_pop;
    N_gen = problem.N_gen;
    t_gen = problem.t_gen;
    obs_Nth = problem.obs_Nth;
    selection_type = problem.selection_type;
      
    %%% GENERATE WRIGHT-FISHER TRAJECTORY IF REQUESTED
    if strcmpi(gen_traj,'wf') || strcmpi(gen_traj,'both')
        
        % Create an observation template
        obs_gen = 1:obs_Nth:N_gen+1;
        N_obs = length(obs_gen);
        % Generate the time points associated with observations
        t_data = t_gen * (0:N_gen);
        t_data = t_data(obs_gen);
        % Prepare storage for Wright-Fisher trajectories
        traj_WF.X = zeros(N_feat,N_obs);
        traj_WF.t = zeros(1,N_obs);
        %  Generate a Wright-Fisher trajectory
        X_trajectory = wrightFisher(N_pop, N_gen, X0, fitness, selection_type);
        % "Observe" this trajectory at observation points and store
        traj_WF.t = t_data;
        traj_WF.X = X_trajectory(:,obs_gen);
        
        % Store in results struct
        trajectories{k}.WF = traj_WF;
        
    end
    
    %%% GENERATE REPLICATOR TRAJECTORY IF REQUESTED
    if strcmpi(gen_traj,'rep') || strcmpi(gen_traj,'both')
        
        % Create a vector of timepoints to use in storing ODE simulation output
        timepoints = linspace(0,N_gen*t_gen,N_simpts);
        % Prepare a mass matrix to confine results to probability simplex
        M = eye(N_feat);
        M(N_feat,N_feat) = 0;
        % Initialise ODE solve tolerance and most strict tolerance
        tol = 1e-3;
        min_tol = 1e-10;
        running = true;
        success = false;
        
        % Disable warnings as we will try many tolerances
        warning('off','MATLAB:ode15s:IntegrationTolNotMet');
        
        % Solve the replicator ODE
        while running

            % Attempt solve using the current tolerance       
            odeoptions = odeset('Mass',M,'AbsTol',tol,'RelTol',tol);
            [T_rep,X_rep] = ode15s( @(t,X) 1/t_gen * replicatorRHS(X,fitness,selection_type), timepoints, X0, odeoptions);
                                   
            % Terminate if integration successful or tolerance too low
            if abs(T_rep(end) - timepoints(end)) < 1e-10 || tol <= min_tol
                running = false;
                success = true;
            end
            
            % Terminate with failure if tolerance too small
            if tol <= min_tol && ~success
                running = false;
                fprintf('\n WARNING: There was a failure to solve the ODE even with strictest tolerance! Output will be truncated.\n');
            end
                
            % Otherwise, decrease tolerance and try again
            tol = tol * 0.01;
            
        end
                  
        % Re-enable warnings to not bother the user's session
        warning('on','MATLAB:ode15s:IntegrationTolNotMet');
        
        % Store the trajectory
        traj_rep.t = T_rep;
        traj_rep.X = X_rep';
        
        % Store in results struct
        trajectories{k}.rep = traj_rep;
        
    end
        
end

% If the input was not a cell array, remove output from cell
if ~cell_input
    trajectories = trajectories{1}; 
end
    
end