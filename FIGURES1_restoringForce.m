function FIGURES1_restoringForce()
% This function demonstrates the restoring force effect that causes MLE
% parameter estimates to be incorrect for the s = 0 case

%%% INITIALISATION

% Load in the basic allele problem
load('problems_basic.mat','allele_basic');
problem = allele_basic;
problem.obs_Nth = 1;
problem.N_pop = 1000;
problem.selection_type = 2;
problem.t_gen = 1;
problem.X0 = [0.2;0.8];
N_realisations = 200;

% Specify the selection strength and dominance
s_true = 0;
h_true = 0.8;

% Specify the number of tiles to use per plot (higher numbers gives a
% smaller gap between the two plots
N_tiles = 15;

% Set random seed - default value is 7 but this one was adjusted to
% generate a nice example where drift in data doesn't look like selection
rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare fitness function and the problem's actual fitness
fitness_fun = @(X,s,hs) [1 + s, 1 + hs; 1 + hs, 1] * X;
problem.fitness = @(X) fitness_fun(X,s_true,h_true*s_true);

% Create data with no selection present
trajWF = simulateTrajectories(problem,'WF');
trajWF = trajWF.WF;

% Fit maximum likelihood to this
optimiser_options = optimset('fminunc');
optimiser_options.MaxFunEvals = 20000;
optimiser_options.MaxIter = 10000;
params = fminunc(@(params) -WFLikelihood2(trajWF, setfield(problem, 'fitness', @(X) fitness_fun(X, params(1), params(2)))), [s_true;h_true*s_true], optimiser_options);

% Prepare figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
tiles_obj = tiledlayout(1,2*(N_tiles+1) - 1);
tiles_obj.TileSpacing = 'compact';
tiles_obj.Padding = 'compact';

% Initialise the plot limits
cur_min = Inf;
cur_max = -Inf;

% Plot many Wright-Fisher trajectories using the true fitness
nexttile([1 N_tiles]);
hold on;
for k = 1:N_realisations
    
    % Generate and plot trajectory
    traj_here = simulateTrajectories(problem,'WF');
    plot(traj_here.WF.t, traj_here.WF.X(1,:), '-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    
    % Also track limits of the trajectories and update if needed
    if min(traj_here.WF.X(1,:)) < cur_min
        cur_min = min(traj_here.WF.X(1,:));
    end
    if max(traj_here.WF.X(1,:)) > cur_max
        cur_max = max(traj_here.WF.X(1,:));
    end
    
end

% Plot the data over the top
plot(trajWF.t, trajWF.X(1,:), 'k.', 'MarkerSize', 30);

% Add axes labels and improve axes appearance
xlabel('Number of Generations','FontSize',24);
ylabel('Allele Frequency, x_1','FontSize',24);
set(gca,'FontSize',24,'LineWidth',2);

% Add a title
title('No Selection','FontSize',30, 'Interpreter', 'LaTeX');

% Plot a blank figure for gap
nexttile;
set(gca,'Visible',false);

% Plot many Wright-Fisher trajectories using the learned fitness
nexttile([1 N_tiles]);
hold on;
for k = 1:N_realisations
    
    % Generate and plot trajectory
    traj_here = simulateTrajectories(setfield(problem, 'fitness', @(X) fitness_fun(X, params(1), params(2))),'WF');
    plot(traj_here.WF.t, traj_here.WF.X(1,:), '-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    
    % Also track limits of the trajectories and update if needed
    if min(traj_here.WF.X(1,:)) < cur_min
        cur_min = min(traj_here.WF.X(1,:));
    end
    if max(traj_here.WF.X(1,:)) > cur_max
        cur_max = max(traj_here.WF.X(1,:));
    end
    
end

% Plot the data over the top
plot(trajWF.t, trajWF.X(1,:), 'k.', 'MarkerSize', 30);

% Add axes labels and improve axes appearance
xlabel('Number of Generations','FontSize',24);
ylabel('Allele Frequency, x_1','FontSize',24);
set(gca,'FontSize',24,'LineWidth',2);

% Add a title
title(['$s$ = ',num2str(params(1)),', $h$ = ',num2str(params(2)/params(1))],'FontSize',30, 'Interpreter', 'LaTeX');

% Set axis limits to be consistent for both
nexttile(1);
ylim([cur_min cur_max]);
nexttile(N_tiles+2);
ylim([cur_min cur_max]);
