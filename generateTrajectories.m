function trajectories = generateTrajectories( problems )
% This function generates both a Wright-Fisher trajectory and the
% replicator trajectory for the problem (or set of problems if given a cell
% array) that are provided as input. Part of these problem objects should
% also specify how many generations to run the dynamics for, and how
% frequently these dynamics (an integer number of generations) should be 
% observed

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
    
    % Create a vector of timepoints to use in storing ODE simulation output
    timepoints = linspace(0,max(t_data),N_simpts);
    % Simulate the given dynamics' replicator equation
    [T_rep,X_rep] = ode15s( @(t,X) 1/t_gen * replicatorRHS(X,fitness,selection_type), timepoints, X0);
    traj_rep.t = T_rep;
    traj_rep.X = X_rep';
    
    % Store this trajecory in the output array
    trajectories{k} = struct('WF',traj_WF,'rep',traj_rep);
    
end

% If the input was not a cell array, remove output from cell
if ~cell_input
    trajectories = trajectories{1}; 
end
    
end