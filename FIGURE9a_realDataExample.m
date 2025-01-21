function FIGURE9a_realDataExample()
% This function simply creates a visual demonstration of the equation
% learning process and its performance

%%% INITIALISATION

% Load in the panaxia data
load('panaxia_data.mat','panaxia_data');
% Store starting year (for shifting data to/from 0)
start_yr = panaxia_data.Year(1);
% Time values are years, shifted to start at zero
data.t = panaxia_data.Year - start_yr;
% X values are the frequency of the unfavoured allele, then favoured
data.X = [panaxia_data.MedionigraAlleleFreq'; panaxia_data.RegularAlleleFreq'];

% Set up the problem object associated with the replicator plotting
problem.N_gen = max(data.t) + 1;
problem.t_gen = 1;
problem.selection_type = 2;
problem.N_feat = 2;
problem.X0 = data.X(:,1);

% Prepare equation learning options
ELoptions.library_orders = 1;
ELoptions.deriv_method = 'gp';
ELoptions.use_smoothed = false;
ELoptions.symmetric_payoff = true;

% Perform equation learning using gradient matching
[library, fitfun] = evolutionLearning(data, problem.selection_type, setfield(ELoptions, 'regression_type', 'grad'));

% Output the parameters estimated
A = libraryToPayoff(library,2);
params = extractBialleleParams(A,problem.selection_type);
grad_s = params.s;
grad_h = params.hs / params.s;

% Simulate the trajectory for the learned fitness
learned_traj = simulateTrajectories( setfield(problem,'fitness',@(X) evaluateLibrary(X,library)), 'rep');

% Also do a separate learning of the actual data directly
library = evolutionLearning(data, problem.selection_type, setfield(ELoptions, 'regression_type', 'ls'));

% Output the parameters estimated
A = libraryToPayoff(library,2);
params = extractBialleleParams(A,problem.selection_type);
ls_s = params.s;
ls_h = params.hs / params.s;

% Simulate the trajectory for the learned fitness
nlin_traj = simulateTrajectories( setfield(problem,'fitness',@(X) evaluateLibrary(X,library)), 'rep');


%%% PLOTTING

% Prepare the figure
figure('units','Normalized','Position',[0.2 0.2 0.2 0.2]);
hold on;

% Plot the data
plot(start_yr + data.t, data.X(1,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);
legend_txt{1} = 'Observed Data';

% Generate the fitted function
Ffit = fitfun{1}.F(data.t);

% Plot the fitted function used for derivative estimates (undo the
% non-dimensionalisation)
plot(start_yr + data.t, Ffit(1,:), 'k--', 'LineWidth', 2, 'color', [0.5 0.5 0.5]);
legend_txt{end+1} = 'Derivative Estimator';

% Plot the learned replicator
indata = learned_traj.rep.t <= max(data.t);
plot(start_yr + learned_traj.rep.t(indata), learned_traj.rep.X(1,indata), 'LineWidth', 3, 'color', [0 0 1]);
legend_txt{end+1} = 'Gradient Matched Replicator';
plot(start_yr + learned_traj.rep.t(~indata), learned_traj.rep.X(1,~indata), '--', 'LineWidth', 3, 'color', [0 0 1]);
legend_txt{end+1} = '';

% Plot the nlinfit replicator
indata = nlin_traj.rep.t <= max(data.t);
plot(start_yr + nlin_traj.rep.t(indata), nlin_traj.rep.X(1,indata), 'LineWidth', 3, 'color', [1 0 0]);
legend_txt{end+1} = 'Least-squares Replicator';
plot(start_yr + nlin_traj.rep.t(~indata), nlin_traj.rep.X(1,~indata), '--', 'LineWidth', 3, 'color', [1 0 0]);
legend_txt{end+1} = '';

% Find the maximum X value across data and curves to be plotted
max_X = max([ max(data.X(1,:)); max(learned_traj.rep.X(1,:)); max(nlin_traj.rep.X(1,:))]);

% Add the curve labels (can be manually dragged into best position if required)
text_t = start_yr + 0.65*max(data.t);
text_X = max_X * 0.7;
text(text_t, text_X+0.025*max_X, ['s = ',num2str(grad_s,'%0.2f'),', h =',num2str(grad_h,'%0.2f')], 'FontSize', 24, 'Color', [0 0 1]);
text(text_t, text_X-0.025*max_X, ['s = ',num2str(ls_s,'%0.2f'),', h = ',num2str(ls_h,'%0.2f')], 'FontSize', 24, 'Color', [1 0 0]);

% Re-plot the data to ensure it appears on top
plot(start_yr + data.t, data.X(1,:), '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0]);

% Add the legend
legend(legend_txt, 'FontSize', 24,'Box','off');

% Clean up the plot
axis( [ [start_yr, start_yr+max(data.t)], [0, 1.05*max_X] ] );
axis square;
xlabel('Year', 'FontSize', 24);
ylabel('Medionigra Allele Frequency', 'FontSize', 24);
set(gca,'LineWidth',2,'FontSize',20);
