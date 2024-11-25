function FIGURE3_basicAlleleExamples()
% This function simply creates a visual demonstration of the equation
% learning process and its performance

%%% INITIALISATION

% Load in the basic allele problem
load('problems_basic.mat','allele_basic');
problem = allele_basic;
problem.obs_Nth = 1;
problem.N_pop = 1000;
problem.selection_type = 2;
problem.t_gen = 1;
obs_gen = 20;
sim_gen = 50;
N_reps = 1;

% Specify the seeds for RNG to plot with (default value of 7 plus one more)
seed_nums = [7,20];


% Generate an example plot for each random seed provided
for s = 1:length(seed_nums)
    
    % Set random seed
    rng( seed_nums(s) );
    
    
    %%% DATA GENERATION AND EQUATION LEARNING
    
    % Generate the Wright-Fisher data for this problem
    problem.N_gen = obs_gen;
    for k = 1:N_reps
        trajWF = generateTrajectories(problem,'WF');
        trajWFs{k} = struct('t',trajWF.WF.t,'X',trajWF.WF.X);
        trajWFs_nonD{k} = struct('t',trajWF.WF.t / problem.t_gen, 'X', trajWF.WF.X);
    end
    
    % Perform equation learning
    ELoptions.library_orders = 0:1;
    ELoptions.deriv_method = 'gp';
    ELoptions.use_smoothed = false;
    ELoptions.symmetric_payoff = false;
    [library, fitfun] = evolutionLearning(trajWFs_nonD, problem.selection_type, setfield(ELoptions, 'regression_type', 'grad'));
    
    % Update the problem length so we can plot over a longer time-frame
    problem.N_gen = sim_gen;
    
    % Generate the replicator data using true fitness over this longer period
    traj = generateTrajectories(problem,'rep');
    
    % Create a copy of the problem but with the learned fitness function
    learned_problem = problem;
    learned_problem.fitness = @(X) evaluateLibrary(X,library);
    learned_traj = generateTrajectories(learned_problem, 'rep');
    
    % Also do a separate learning of the actual data directly
    library = evolutionLearning(trajWFs_nonD, problem.selection_type, setfield(ELoptions, 'regression_type', 'ls'));
    nlin_problem = problem;
    nlin_problem.fitness = @(X) evaluateLibrary(X,library);
    nlin_traj = generateTrajectories(nlin_problem, 'rep');
    
    % Evaluate the fitted function used in gradient matching for plotting
    Ffit = fitfun{1}.F(traj.rep.t(traj.rep.t <= obs_gen*problem.t_gen) / problem.t_gen);
    
    
    %%% PLOTTING
    
    % Prepare the figure
    figure('units','Normalized','Position',[0.2 0.2 0.2 0.2]);
    hold on;
    
    % Initialise legend text
    legend_txt = cell(0);
    
    % Plot the data
    for k = 1:N_reps
        plot(trajWFs{k}.t, trajWFs{k}.X(1,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);
        if k == 1
            legend_txt{k} = 'Data';
        else
            legend_txt{k} = '';
        end
    end
    
    % Plot the true replicator
    plot(traj.rep.t, traj.rep.X(1,:), 'LineWidth', 4, 'color', [0.4 0.4 0.4]);
    legend_txt{end+1} = 'True Selective Pressure';
    
    % Plot the least squares replicator
    plot(nlin_traj.rep.t(traj.rep.t <= obs_gen*problem.t_gen), nlin_traj.rep.X(1,traj.rep.t <= obs_gen*problem.t_gen), 'LineWidth', 3, 'color', [0.8 0 0]);
    legend_txt{end+1} = 'Nonlinear Least Squares';
    plot(nlin_traj.rep.t(traj.rep.t > obs_gen*problem.t_gen), nlin_traj.rep.X(1,traj.rep.t > obs_gen*problem.t_gen), '--', 'LineWidth', 3, 'color', [0.8 0 0]);
    legend_txt{end+1} = '';
    
    % Plot the gradient matched replicator
    plot(learned_traj.rep.t(traj.rep.t <= obs_gen*problem.t_gen), learned_traj.rep.X(1,traj.rep.t <= obs_gen*problem.t_gen), 'LineWidth', 3, 'color', [0 0 0.8]);
    legend_txt{end+1} = 'Gradient Matching';
    plot(learned_traj.rep.t(traj.rep.t > obs_gen*problem.t_gen), learned_traj.rep.X(1,traj.rep.t > obs_gen*problem.t_gen), '--', 'LineWidth', 3, 'color', [0 0 0.8]);
    legend_txt{end+1} = '';
    
    % Plot the fitted function used for derivative estimates (undo the
    % non-dimensionalisation)
    plot(traj.rep.t(traj.rep.t <= obs_gen*problem.t_gen), Ffit(1,:), 'k--', 'LineWidth', 3, 'color', [0.7 0.7 0.7]);
    legend_txt{end+1} = 'Derivative Estimator';
    
    
    % Re-plot the data to ensure it appears on top
    for k = 1:N_reps
        plot(trajWFs{k}.t, trajWFs{k}.X(1,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);
    end
    
    % Plot clean-up
    set(gca,'LineWidth',2,'FontSize',20);
    axis square;
    xlabel('Number of Generations','FontSize',20);
    ylabel('Allele Frequency, x_1','FontSize',20);
    % Manually set axis based on known properties of problem
    ylim([0.2 1]);
    
    % Add legend to plot if this is first seed
    if s == 1
        legend(legend_txt, 'FontSize', 24, 'box', 'off');
    end
    
end