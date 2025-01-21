function FIGURES6_RPSManyExamples
%
% This function visualises example trajectories, and the corresponding
% fitted fitness libraries, for the three different rock-paper-scissors 
% scenarios that we use as test problems. Plots are generated at a low
% population level, to show the effects of drift/stochasticity, and also at
% a high population level, to highlight learning of the "true" trajectory


%%% PROBLEM SPECIFICATIONS

% Observation frequency
obs_Nth = 1;
% Number of generations to actually observe
obs_gen = 100;
% Number of generations to simulate trajectories over
sim_gen = 200;
% Number of separate datasets to use in learning
N_data = 1;
% Number of replicates to plot on each triangle
N_reps = 10;
% Population size
N_pops = [1000];
% Selection type
selection_type = 1;

%%% EVOLUTION LEARNING SPECIFICATIONS

ELoptions.library_orders = 0:1;
ELoptions.symmetric_payoff = false;
ELoptions.nlin_silent = true;

%%% INITIALISATION

% Load in the set of tri-allelic selection problems, store them in array
load('problems_RPS.mat','RPS_balanced', 'RPS_attract', 'RPS_repel');
problems_list = {RPS_balanced, RPS_attract, RPS_repel};
problem_names = {'Balanced', 'Attracting', 'Repelling'};


%%% FIGURE CREATION

% Loop over each population size given
for n = 1:length(N_pops)
    
    % Read out this current population size
    N_pop = N_pops(n);
    
    % Specify a random seed for consistency
    rng(7);
    
    % Prepare this figure
    figure('units','Normalized','Position',[0.0 0.0 1.0 1.0]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    
    % Loop over each problem in the set of problems
    for p = 1:length(problems_list)
        
        % Advance tile for plotting of this problem
        nexttile
        
        %%% DATA PREPARATION
        
        % Load in the current problem and apply the above properties
        problem = problems_list{p};
        problem.N_pop = N_pop;
        problem.selection_type = selection_type;
        problem.obs_Nth = obs_Nth;
        
        % Loop over requested number of replicates
        for r = 1:N_reps
            
            % Generate Wright-Fisher data for this problem, including
            % non-dimensionalised form for learning
            problem.N_gen = obs_gen;
            for k = 1:N_data
                trajWF = simulateTrajectories(problem,'WF');
                trajWFs{k} = struct('t',trajWF.WF.t,'X',trajWF.WF.X);
                trajWFs_nonD{k} = struct('t',trajWF.WF.t / problem.t_gen, 'X', trajWF.WF.X);
            end
            
            
            %%% EVOLUTION LEARNING AND TRAJECTORY GENERATION
            
            % Switch to the simulation version of the problem
            problem.N_gen = sim_gen;
            
            % Apply evolution learning to this data, using both types of regression
            library_grad = evolutionLearning(trajWFs_nonD, problem.selection_type, setfield(ELoptions, 'regression_type','grad'));
            library_ls = evolutionLearning(trajWFs_nonD, problem.selection_type, setfield(ELoptions, 'regression_type','ls'));
            
            % Simulate the replicator for the learned fitness functions
            traj_grad = simulateTrajectories(setfield(problem, 'fitness', @(X) evaluateLibrary(X,library_grad)), 'rep' );
            traj_ls = simulateTrajectories(setfield(problem, 'fitness', @(X) evaluateLibrary(X,library_ls)), 'rep' );
            
            % Also simulate the true selective pressure for this problem
            traj_true = simulateTrajectories(problem, 'rep');
            
            
            %%% PLOT ALL CURVES TOGETHER
            
            % Prepare the figure
            hold on;
            
            % Clear legend for each go through loop
            legend_txt = cell(0);
                       
            % Plot the true replicator
            trianglePlot(traj_true.rep.t, traj_true.rep.X', 'LineWidth', 5, 'color', [0.4 0.4 0.4]);
            % Extra legend entries for the triangle
            legend_txt(end+1:end+6) = {''};
            % Add a legend entry for it
            legend_txt{end+1} = 'True Selective Pressure      $\hphantom{f}$';
            
            % Plot the least squares replicator - observation and forecast period
            trianglePlot(traj_ls.rep.t(traj_ls.rep.t <= problem.t_gen * obs_gen), traj_ls.rep.X(:,traj_ls.rep.t <= problem.t_gen * obs_gen)', 'LineWidth', 5, 'color', [1 0.4 0.4]);
            %trianglePlot(traj_ls.rep.t(traj_ls.rep.t > problem.t_gen * obs_gen), traj_ls.rep.X(:,traj_ls.rep.t > problem.t_gen * obs_gen)', ':', 'LineWidth', 5, 'color', [1 0.4 0.4]);
            % Legend entry for curve and blank for forecast
            legend_txt{end+1} = 'Least Squares      $\hphantom{f}$';
            %legend_txt{end+1} = '';
            
            % Plot the gradient matching replicator - observation and forecast period
            trianglePlot(traj_grad.rep.t(traj_grad.rep.t <= problem.t_gen * obs_gen), traj_grad.rep.X(:,traj_grad.rep.t <= problem.t_gen * obs_gen)', 'LineWidth', 5, 'color', [0.4 0.4 1]);
            %trianglePlot(traj_grad.rep.t(traj_grad.rep.t > problem.t_gen * obs_gen), traj_grad.rep.X(:,traj_grad.rep.t > problem.t_gen * obs_gen)', ':', 'LineWidth', 5, 'color', [0.4 0.4 1]);
            % Legend entry for curve and blank for forecast
            legend_txt{end+1} = 'Gradient Matching      $\hphantom{f}$';
            %legend_txt{end+1} = '';
            
            % Re-plot the selective trajectory to ensure it appears on top
            trianglePlot(traj_true.rep.t, traj_true.rep.X', 'LineWidth', 5, 'color', [0.4 0.4 0.4]);
            
            % Add a title to the current plot if this is the top row
            if n == 1
                title(problem_names{p}, 'FontSize', 32, 'Interpreter', 'LaTeX');
            end
            
            % Only actually add the legend to one plot from the last set of
            % figures
            if r == 1 && p == 1 && n == length(N_pops)
                legend(legend_txt, 'FontSize', 28, 'box', 'off', 'Orientation', 'Horizontal', 'Interpreter', 'LaTeX');
            end
            
        end
        
    end
end