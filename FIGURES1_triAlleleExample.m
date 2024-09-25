function FIGURES1_triAlleleExample()
% This function simply creates a visual demonstration of the equation
% learning process and its performance

%%% INITIALISATION

% Load in the basic allele problem
load('problems_inheritance.mat','allele_standard');
problem = allele_standard;
problem.obs_Nth = 25;
problem.N_pop = 5000;
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
ELoptions.deriv_method = 'poly';
ELoptions.kernel_bandwidth = 100;
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

for n = 1:length(problem.X0);

% Prepare the figure
figure('units','Normalized','Position',[0.0 0.0 1.0 1.0]);
hold on;

% Clear legend for each go through loop
legend_txt = cell(0);

% Plot the data
for k = 1:N_reps
    plot(trajWFs{k}.t, trajWFs{k}.X(n,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);
    if k == 1
        legend_txt{k} = 'Data';
    else
        legend_txt{k} = '';
    end
end

% Plot the true replicator
plot(traj.rep.t, traj.rep.X(n,:), 'LineWidth', 4, 'color', [0 0 0]);
legend_txt{end+1} = 'True Selective Pressure';

% Plot the fitted function used for derivative estimates (undo the
% non-dimensionalisation)
plot(traj.rep.t, fitfun.f_funs{n}(traj.rep.t / problem.t_gen), 'k--', 'LineWidth', 2, 'color', [0.5 0.5 0.5]);
legend_txt{end+1} = 'Derivative Estimator';

% Plot the learned replicator
plot(learned_traj.rep.t, learned_traj.rep.X(n,:), 'LineWidth', 3, 'color', [0 0 1]);
legend_txt{end+1} = 'Learned Replicator';

% Plot the nlinfit replicator
plot(nlin_traj.rep.t, nlin_traj.rep.X(n,:), 'LineWidth', 3, 'color', [1 0 0]);
legend_txt{end+1} = 'Least-squares Replicator';

% Re-plot the data to ensure it appears on top
for k = 1:N_reps
    plot(trajWFs{k}.t, trajWFs{k}.X(n,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);
end

legend(legend_txt, 'FontSize', 24);

end