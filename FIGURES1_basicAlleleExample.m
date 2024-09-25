function FIGURES1_basicAlleleExample()
% This function simply creates a visual demonstration of the equation
% learning process and its performance

%%% INITIALISATION

% Load in the basic allele problem
load('problems_basic.mat','allele_basic');
problem = allele_basic;
problem.obs_Nth = 10;
problem.N_pop = 2000;
problem.selection_type = 1;
N_reps = 1;

% Specify a random seed for consistency
rng(7);
 

%%% DATA GENERATION AND EQUATION LEARNING

% Generate the Wright-Fisher and replicator data for this problem
traj = generateTrajectories(problem,'rep');
for k = 1:N_reps
    trajWF = generateTrajectories(problem,'WF');
    trajWFs{k} = struct('t',trajWF.WF.t,'X',trajWF.WF.X);
    trajWFs_nonD{k} = struct('t',trajWF.WF.t / problem.t_gen, 'X', trajWF.WF.X);
end

% Perform equation learning
ELoptions.library_orders = 1;
ELoptions.deriv_method = 'gp';
[fitfun, Flib, K, texts] = evolutionLearning(trajWFs_nonD, problem.selection_type, ELoptions);

% Create a copy of the problem but with the learned fitness function
learned_problem = problem;
learned_problem.fitness = @(X) libraryRHS(X,Flib,K);
learned_traj = generateTrajectories(learned_problem, 'rep');

% Also do a separate learning of the actual data directly
[Flib, K, texts] = LSevolutionLearning(trajWFs_nonD, problem.selection_type, ELoptions);
nlin_problem = problem;
nlin_problem.fitness = @(X) libraryRHS(X,Flib,K);
nlin_traj = generateTrajectories(nlin_problem, 'rep');


%%% PLOTTING

% Prepare the figure
figure('units','Normalized','Position',[0.2 0.2 0.2 0.2]);
hold on;

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
plot(traj.rep.t, traj.rep.X(1,:), 'LineWidth', 4, 'color', [0 0 0]);
legend_txt{end+1} = 'True Selective Pressure';

% Generate the fitted function
Ffit = fitfun{1}.F(traj.rep.t / problem.t_gen);

% Plot the fitted function used for derivative estimates (undo the
% non-dimensionalisation)
plot(traj.rep.t, Ffit(1,:), 'k--', 'LineWidth', 2, 'color', [0.5 0.5 0.5]);
legend_txt{end+1} = 'Derivative Estimator';

% Plot the learned replicator
plot(learned_traj.rep.t, learned_traj.rep.X(1,:), 'LineWidth', 3, 'color', [0 0 1]);
legend_txt{end+1} = 'Learned Replicator';

% Plot the nlinfit replicator
plot(nlin_traj.rep.t, nlin_traj.rep.X(1,:), 'LineWidth', 3, 'color', [1 0 0]);
legend_txt{end+1} = 'Least-squares Replicator';

% Re-plot the data to ensure it appears on top
for k = 1:N_reps
    plot(trajWFs{k}.t, trajWFs{k}.X(1,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);
end

legend(legend_txt, 'FontSize', 24);
