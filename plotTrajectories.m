function plotTrajectories(simResults, problems, pop_level, show_types)
%
% This function plots the trajectories for type I and type II selection, as
% generated by the learned fitnesses
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify how many of the replicates to plot
N_plot = 50;

% Specify the plotting colours
plot_clrs = [ 0.6, 0.0, 0.0;
              1.0, 0.4, 0.4;
              0.0, 0.0, 0.6;
              0.4, 0.4, 1.0  ];

% Loop over type 1 and type 2 selection
for type = show_types
    
    % Extract problem, including its initial condition and true fitness
    problem = problems{type};
    
    % Create also a flipped initial condition problem
    ICproblem = problem;
    ICproblem.X0 = flipud(ICproblem.X0);
    
    % Extract the fitness functions for the two methods for this type
    F_LS = squeeze(simResults.f_mat(pop_level,type,:));
    F_EL = squeeze(simResults.f_mat(pop_level,type+2,:));
    
    % Initialise figure
    figure('units','normalized','Position',[0.2 0.2 0.6 0.6]);
    hold on;
    
    % Loop over the replicates and plot their trajectories
    for k = 1:N_plot
       
        % Create least squares problem 
        LS_problem = problem;
        LS_problem.fitness = F_LS{k};
        
        % Generate and plot this trajectory
        traj = generateTrajectories(LS_problem,'rep');
        plot(traj.rep.t, traj.rep.X(1,:), 'LineWidth', 2, 'Color', plot_clrs(1,:));
        
        % Create flipped least squares problem 
        LS_problem = ICproblem;
        LS_problem.fitness = F_LS{k};
        
        % Generate and plot this trajectory
        traj = generateTrajectories(LS_problem,'rep');
        plot(traj.rep.t, traj.rep.X(1,:), 'LineWidth', 2, 'Color', plot_clrs(2,:));
        
        % Create evolution learning problem 
        EL_problem = problem;
        EL_problem.fitness = F_EL{k};
        
        % Generate and plot this trajectory
        traj = generateTrajectories(EL_problem,'rep');
        plot(traj.rep.t, traj.rep.X(1,:), 'LineWidth', 2, 'Color', plot_clrs(3,:));
        
        % Create flipped evolution learning problem 
        EL_problem = ICproblem;
        EL_problem.fitness = F_EL{k};
        
        % Generate and plot this trajectory
        traj = generateTrajectories(EL_problem,'rep');
        plot(traj.rep.t, traj.rep.X(1,:), 'LineWidth', 2, 'Color', plot_clrs(4,:));
        
    end
    
    % Generate and plot the true trajectories
    traj = generateTrajectories(problem,'rep');
    plot(traj.rep.t, traj.rep.X(1,:), 'k', 'LineWidth', 3);
    traj = generateTrajectories(ICproblem,'rep');
    plot(traj.rep.t, traj.rep.X(1,:), 'k', 'LineWidth', 3);

end