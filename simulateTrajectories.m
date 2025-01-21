function trajectories = simulateTrajectories( problems, trajectory_type )
%
%     trajectories = simulateTrajectories(problem)
%     trajectories = simulateTrajectories(problem, trajectory_type)
%     trajectories = simulateTrajectories({problem1,problem2,...})
%     trajectories = simulateTrajectories({problem1,problem2,...}, trajectory_type)
%
% This function generates a Wright-Fisher and/or a replicator trajectory
% for the input problem object, or list of problem objects. If the user
% does not specify which type of trajectory to generate, both will be
% produced.
%
% Problem objects not only specify the selection dynamics that are active,
% but also the number of generations to simulate over, the length of time
% each generation physically represents, frequency with which observations
% should be recorded, and so on.
%
%  INPUTS:
%
%         problems - A problem object, which specifies the fitness function
%                    and type of selection to apply, as well as observation
%                    period, observation frequency, etc.
%                    Example problem objects can be seen in the function
%                    createTestProblems.m
%
%  trajectory_type - The type of trajectory or trajectories to generate.
%                    Valid options include
%                        'wf': Wright-Fisher trajectory over the specified
%                              number of generations, with each N_obs-th
%                              generation being recorded, and using a 
%                              population size of N_pop
%                       'rep': Replicator trajectory for this problem(s),
%                              over the same period of time for which the
%                              Wright-Fisher would be generated. The
%                              observations are generated at the specified
%                              number of points, below
%                      'both': Simulates and stores both 'wf' and 'rep'
%                              
%
%  OUTPUTS:
%
%     trajectories -  A structure that contains the (t,X) data for the
%                     trajectory or trajectories requested. These are
%                     of form  (trajectories.WF.t, trajectories.WF.X)
%                      and/or  (trajectories.rep.t, trajectories.rep.X)
%                     If a list of problem objects was provided, the output
%                     will be a list of trajectory objects, each as just
%                     described.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PARAMETER DEFINITIONS

% Number of points to use for replicator trajectories
N_simpts = 10001;
% Initial tolerance to trial replicator solves
init_tol = 1e-7;
% Minimum tolerance for replicator solves before abandoning
min_tol = 1e-11;
% Tolerance decrease factor (scales tolerance upon integration failure)
tol_factor = 0.01;
% Flag for solving using either a standard or DAE formulation
solve_as_DAE = true;
 
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
    trajectory_type = 'both';
end

% Loop over each problem, creating the trajectory object for each
trajectories = cell(N_prob,1);
for k = 1:N_prob
   
    % Read out all the current problem's settings
    problem = problems{k};
    N_feat = problem.N_feat;
    fitness = problem.fitness;
    X0 = problem.X0;
    N_gen = problem.N_gen;
    t_gen = problem.t_gen;
    selection_type = problem.selection_type;
      
    %%% GENERATE WRIGHT-FISHER TRAJECTORY IF REQUESTED
    if strcmpi(trajectory_type,'wf') || strcmpi(trajectory_type,'both')
        
        % Extract additional WF-specific information
        N_pop = problem.N_pop;
        obs_Nth = problem.obs_Nth;
        
        % Create an observation template
        obs_gen = 1:obs_Nth:N_gen+1;
        N_obs = length(obs_gen);
        % Generate the time points associated with observations
        t_data = t_gen * (0:N_gen);
        t_data = t_data(obs_gen);
        % Prepare storage for Wright-Fisher trajectories
        traj_WF.t = zeros(1,N_obs);
        traj_WF.X = zeros(N_feat,N_obs);
        % Generate a Wright-Fisher trajectory
        X_trajectory = wrightFisher(N_pop, N_gen, X0, fitness, selection_type);
        % "Observe" this trajectory at observation points and store
        traj_WF.t = t_data;
        traj_WF.X = X_trajectory(:,obs_gen);
        
        % Store in results struct
        trajectories{k}.WF = traj_WF;
        
    end
    
    %%% GENERATE REPLICATOR TRAJECTORY IF REQUESTED
    if strcmpi(trajectory_type,'rep') || strcmpi(trajectory_type,'both')
        
        % Create a vector of timepoints to use in storing ODE simulation output
        timepoints = linspace(0,N_gen*t_gen,N_simpts);
        
        % Create a mass matrix if we will be using one
        if solve_as_DAE
            M = diag([ones(N_feat-1,1);0]);
        end
        
        % Disable warnings as we will try many tolerances
        warning('off','MATLAB:ode15s:IntegrationTolNotMet');
        
        % Initialise tolerance
        tol = init_tol;
        
        % Constantly attempt solves until success, or tolerance too small
        running = true;
        success = false;
        while running
            
            % Attempt a DAE solve if the flag to do so is set
            DAE_failed = false;
            if solve_as_DAE
                
                % Attempt solve at current tolerance
                try
                    odeoptions = odeset('Mass',M,'AbsTol',tol,'RelTol',tol);
                    [T_rep, X_rep] = ode15s( @(t,X) 1/t_gen * replicatorRHS(X,fitness,selection_type,true), timepoints, X0, odeoptions );
                catch
                    DAE_failed = true;
                end
                
                % If simulation did not make it to final time, mark failure
                if ~DAE_failed && abs(T_rep(end) - timepoints(end)) > 1e-10
                    DAE_failed = true;
                end
                
            end
            
            % If DAE solve failed, or not trying, solve using standard ODE
            if ~solve_as_DAE || DAE_failed
                odeoptions = odeset('AbsTol',tol,'RelTol',tol);
                [T_rep, X_rep] = ode15s( @(t,X) 1/t_gen * replicatorRHS(X,fitness,selection_type,false), timepoints, X0, odeoptions);
            end
            
            % Terminate if integration successful
            if abs(T_rep(end) - timepoints(end)) <= 1e-10
                running = false;
                success = true;
            end
            
            % Terminate without success if tolerance too small
            if tol <= min_tol && ~success
                running = false;
                fprintf('\n WARNING: There was a failure to solve the ODE even with strictest tolerance! Output will be truncated.\n');
            end
        
            % Otherwise, decrease tolerance and try again
            tol = tol * tol_factor;
            
        end
                
        % Store the trajectory
        traj_rep.t = T_rep;
        traj_rep.X = X_rep';
        
        % Store in results struct
        trajectories{k}.rep = traj_rep;
        
        % Re-enable warnings to not bother the user's session
        warning('on','MATLAB:ode15s:IntegrationTolNotMet');
        
    end
        
end

% If the input was not a cell array, remove output from cell
if ~cell_input
    trajectories = trajectories{1}; 
end
    
end