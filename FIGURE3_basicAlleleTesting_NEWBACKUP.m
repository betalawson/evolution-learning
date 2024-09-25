function FIGURE3_basicAlleleTesting_NEWBACKUP(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying selection strength and
% dominance for a bi-allelic inheritance scenario

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify whether to save figures or just leave them all displayed
save_figs = false;

% Specify the population level to compare methods for
pop_level = 2;

test_types = 1:2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INPUT HANDLING

% Assume not regenerating data unless specifically requested
if nargin < 1
    regenerate = false;
elseif regenerate
    userYN = input("All data will be regenerated. Please confirm Y/N: ",'s');
    if strcmpi(userYN,'y')
        regenerate = true;
    else
        regenerate = false;
    end
end


%%% PROBLEM PREPARATION

% Load in the basic problem - type I and type II selection versions.
% Settings for these will be adjusted as needed
load('problems_basic.mat','allele_basic');

% Create a type I and type II selection version of the problem
problems = cell(1,2);
for type = 1:2
    problem = allele_basic;
    problem.selection_type = type;
    problems{type} = problem;
end

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 5;

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2], [1] };
symmetric_flags = { false, false, false, true };

% Set up the filenames where data will be saved/loaded from
data_prefix = 'DATA_BIALL_orders';
data_methodnames = {'orders01','orders1','orders012','symm'};

% Count provided number of methods
N_methods = length(data_methodnames);


%%% GENERATE DATA IF IT IS NOT PRESENT (OR IF ASKED)

% General linear equation learning
title_txt = '$p \leq 1$';
data_filename = 'DATA_Npop_orders01';
% Regenerate if requested or data not present
if ~exist([data_filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.library_orders = [0,1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
else
    load([data_filename,'.mat'],'simResults');
end

% This is the first problem so use it to set up the population counts
N_pops = length(simResults.pop_sizes);
x_txts = cell(1,N_pops);
for k = 1:N_pops
    x_txts{k} = ['10^{',num2str(log10(simResults.pop_sizes(k))),'}'];
end
x_txts{N_pops+1} = '\infty';

% Actually plot
fig_filenames = {'FIGS1a_biallele_pops-lt1_s', 'FIGS1b_biallele_pops-lt1_h', 'FIGS1c_biallele_pops-lt1_l2'};
plotBoxPlot(simResults, x_txts, title_txt, fig_filenames, save_figs, true);

% Strictly payoff matrix elements
title_txt = '$p = 1$';
data_filename = 'DATA_Npop_orders1';
% Regenerate if requested or data not present
if ~exist([data_filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.library_orders = [1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
else
    load([data_filename,'.mat'],'simResults');
end    
fig_filenames = {'FIG3a_biallele_pops_s', 'FIG3b_biallele_pops_h', 'FIG3c_biallele_pops_l2'};
plotBoxPlot(simResults, x_txts, title_txt, fig_filenames, save_figs, true);

% Quadratic fitness function
title_txt = '$p \leq 2$';
data_filename = 'DATA_Npop_orders012';
% Regenerate if requested or data not present
if ~exist([data_filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.library_orders = [0,1,2];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
else
    load([data_filename,'.mat'],'simResults');
end    
fig_filenames = {'FIGS1d_biallele_pops-lt2_s', 'FIGS1e_biallele_pops-lt2_h', 'FIGS1f_biallele_pops-lt2_l2'};
plotBoxPlot(simResults, x_txts, title_txt, fig_filenames, save_figs, true);

% Specifically a symmetric payoff matrix
title_txt = 'Symm.';
data_filename = 'DATA_Npop_symm';
% Regenerate if requested or data not present
if ~exist([data_filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = true;
    ELoptions.library_orders = [1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
else
    load([data_filename,'.mat'],'simResults');
end    
fig_filenames = {'FIGS1g_biallele_pops-symm_s', 'FIGS1h_biallele_pops-symm_h', 'FIGS1i_biallele_pops-symm_l2'};
plotBoxPlot(simResults, x_txts, title_txt, fig_filenames, save_figs, true);


%%% EFFECTS OF METHOD CHOICE (FOR A SINGLE POPULATION SIZE)

% Define the list of filenames to extract data from
data_filenames = {'DATA_Npop_orders1', 'DATA_Npop_orders01', 'DATA_Npop_orders012', 'DATA_Npop_symm'};

% Define the methods text used on the figure
method_txts = {'p = 1', 'p \leq 1', 'p \leq 2', 'Symm.'};

% Load in the first dataset to set up the structure
load([data_filenames{1},'.mat'],'simResults');
s_mat = zeros([length(data_filenames), 4, size(simResults.s_mat,3)]);
h_mat = zeros([length(data_filenames), 4, size(simResults.s_mat,3)]);
l2_mat = zeros([length(data_filenames), 4, size(simResults.s_mat,3)]);
s_mat(1,:,:) = simResults.s_mat(pop_level,:,:);
h_mat(1,:,:) = simResults.h_mat(pop_level,:,:);
l2_mat(1,:,:) = simResults.l2_mat(pop_level,:,:);

% Load in each other dataset and extract results for this population level
for m = 2:length(data_filenames)
    load([data_filenames{m},'.mat'],'simResults');
    s_mat(m, :, :) = simResults.s_mat(pop_level,:,:);
    h_mat(m, :, :) = simResults.h_mat(pop_level,:,:);
    l2_mat(m, :, :) = simResults.l2_mat(pop_level,:,:);
end

% Plot these results
fig_filenames = {'FIG3d_biallele_methods_s', 'FIG3e_biallele_methods_h', 'FIG3f_biallele_methods_l2'};
compiledResults = struct('s_mat',s_mat,'h_mat',h_mat,'l2_mat',l2_mat,'s_true',simResults.s_true,'h_true',simResults.h_true);
plotBoxPlot(compiledResults, method_txts, ['$N = ',num2str(simResults.pop_sizes(pop_level)),'$'], fig_filenames, save_figs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, filename)
%
% This subfunction generates data for equation learning applied to the
% provided allele inheritance problems, returning the estimated values for
% the parameters (s,h) for different population sizes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify a random seed for consistency - this should give the same
% trajectories for each different approach to match, also
rng(7);

% Loop over the list of population sizes, learning parameters s and h for
% each
N_pops = length(pop_sizes);
fLS_mat = cell(N_pops, 2, N_rep);
fEL_mat = cell(N_pops, 2, N_rep);
sLS_mat = zeros(N_pops, 2, N_rep);
hLS_mat = zeros(N_pops, 2, N_rep);
l2LS_mat = zeros(N_pops, 2, N_rep);
sEL_mat = zeros(N_pops, 2, N_rep);
hEL_mat = zeros(N_pops, 2, N_rep);
l2EL_mat = zeros(N_pops, 2, N_rep);
timing = struct('EL',0,'LS',0);
for k = 1:N_pops
   
    % Loop over the two types of selection and generate results for each
    for type = test_types
    
        % Overwrite the population size for this problem
        problem = problems{type};
        problem.N_pop = pop_sizes(k);
    
        % Initialise timing storage
        LS_toc = zeros(1,N_rep);
        EL_toc = zeros(1,N_rep);
        
        % Repeatedly generate new Wright-Fisher data and learn s and h
        parfor r = 1:N_rep
        
            % Generate the Wright-Fisher and replicator data
            traj = generateTrajectories(problem);
        
            % Convert to non-dimensionalised time for equation learning
            WF_nonD = struct('t',traj.WF.t / problem.t_gen, 'X', traj.WF.X);

            % Perform equation learning
            tic;
            [Flib, K, texts] = LSevolutionLearning(WF_nonD, type, ELoptions);
            LS_toc(r) = toc;
            
            % Store the fitness function in matrix
            fLS_mat{k,type,r} = @(X) libraryRHS(X, Flib, K);
            
            % Extract payoff matrix and normalise using invariant operations
            A = fitnessToPayoff(K, texts, 2);
            [s,h] = extract2AlleleParams(A, type);

            % Extract the values of s and h and store them 
            sLS_mat(k,type,r) = s;
            hLS_mat(k,type,r) = h;
            
            % Simulate the learned dynamics and evaluate the L2 error
            % between the learned dynamics and the actual replicator
            learned_problem = problem;
            learned_problem.fitness = @(X) libraryRHS(X,Flib,K);
            learned_traj = generateTrajectories(learned_problem,'rep');
            l2LS_mat(k,type,r) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));   
                        
            % Now use our gradient matching method
            tic;
            [fitfun, Flib, K, texts] = evolutionLearning(WF_nonD, type, ELoptions);
            EL_toc(r) = toc;
            
            % Store the fitness function in matrix
            fEL_mat{k,type,r} = @(X) libraryRHS(X, Flib, K);
            
            % Extract payoff matrix and normalise using invariant operations
            A = fitnessToPayoff(K, texts, 2);
            [s,h] = extract2AlleleParams(A, type);

            % Extract the values of s and h and store them 
            sEL_mat(k,type,r) = s;
            hEL_mat(k,type,r) = h;
            
            % Simulate the learned dynamics and evaluate the L2 error
            % between the learned dynamics and the actual replicator
            learned_problem = problem;
            learned_problem.fitness = @(X) libraryRHS(X,Flib,K);
            learned_traj = generateTrajectories(learned_problem,'rep');
            l2EL_mat(k,type,r) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));                       
             
        end
        
        timing.LS = timing.LS + sum(LS_toc);
        timing.EL = timing.EL + sum(EL_toc);
        
    end    
    
end

% Now gather together the data
s_mat(:,1:2,:) = sLS_mat;
s_mat(:,3:4,:) = sEL_mat;
h_mat(:,1:2,:) = hLS_mat;
h_mat(:,3:4,:) = hEL_mat;
l2_mat(:,1:2,:) = l2LS_mat;
l2_mat(:,3:4,:) = l2EL_mat;
f_mat(:,1:2,:) = fLS_mat;
f_mat(:,3:4,:) = fEL_mat;


% Now generate "perfect" replicator data with known derivative and estimate
% using it for reference
s_perfect = zeros(1,4);
h_perfect = zeros(1,4);
l2_perfect = zeros(1,4);
for type = test_types
    
    % Load in problem
    problem = problems{type};
    
    % Generate the replicator data, including derivatives
    traj = generateTrajectories(problem);
    X_true = traj.rep.X;
    Xdash_true = zeros(size(X_true));
    for j = 1:size(X_true,2)
        Xdash_true(:,j) = replicatorRHS(X_true(:,j),problem.fitness,type);
    end
    
    % Calculate the best-fit parameters 
    [Flib, K, texts] = LSevolutionLearning(struct('t',traj.rep.t/problem.t_gen,'X',traj.rep.X), type, ELoptions);
    
    % Extract payoff matrix and normalise using invariant operations
    A = fitnessToPayoff(K,texts,2);
    [s_perfect(type),h_perfect(type)] = extract2AlleleParams(A,type);
    
    % Find the L2 error to the true selective trajectory
    learned_problem = problem;
    learned_problem.fitness = @(X) libraryRHS(X,Flib, K);
    learned_traj = generateTrajectories(learned_problem);
    l2_perfect(type) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));
    
    % Provide the perfect data and use it for evolution learning
    ELoptions.Xdash_data = Xdash_true;
    [fitfun, Flib, K, texts] = evolutionLearning(traj.rep, type, ELoptions);
    
    % Extract payoff matrix and normalise using invariant operations
    A = fitnessToPayoff(K, texts, 2);
    [s_perfect(type+2),h_perfect(type+2)] = extract2AlleleParams(A, type);

    % Find the L2 error to the true selective trajectory
    learned_problem = problem;
    learned_problem.fitness = @(X) libraryRHS(X,Flib, K);
    learned_traj = generateTrajectories(learned_problem);
    l2_perfect(type+2) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));
    
end

% Extract the true values of s and h from the problem
A = problems{1}.fitness(eye(length(problem.X0)));    
[s_true, h_true] = extract2AlleleParams(A,1);

% Store all data in a single results object and save it
simResults = struct('f_mat',{f_mat},'s_mat',s_mat,'h_mat',h_mat,'l2_mat',l2_mat,'s_perfect',s_perfect,'h_perfect',h_perfect,'l2_perfect',l2_perfect,'s_true',s_true,'h_true',h_true,'pop_sizes',pop_sizes,'timing',timing);
save([filename,'.mat'],'simResults');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function plotBoxPlot(simResults, x_txts, title_txt, fig_filenames, save_figs, plot_pop)
%
% This subfunction plots the results of equation learning on a boxplot,
% showing the results across repeated replicates as contained in the
% provided structure, 'simResults'. 'x_txts' and 'title_txt' define the
% text to use on x-axis ticks, and additional title text, respectively.
% The figures will be saved as .eps files, unless the flag below is turned
% off.
%
% If the optional argument 'plot_pop' is provided as true, then the x-axis
% is given a label 'N', and the results for perfect data are also shown
% (provide the label for this in 'x_txts')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the colours to use in plotting the two selection types
plot_clrs = [   [ 0.6, 0.0, 0.0 ];          % Dark red
                [ 0.0, 0.0, 0.6 ];          % Dark blue
                [ 1.0, 0.4, 0.4 ];          % Light red
                [ 0.4, 0.4, 1.0 ]  ];       % Light blue
          
% Specify the fraction of data to keep
keep_frac = 0.95;

% Specify the minimum and maximum allowable values for s and h on the plots
% These are used where the range of 'keep_frac' proportion of the data
% becomes too large
y_absmin = -0.75;
y_absmax = 1.75;

% Figure dimensions
fig_x = 0.3;
fig_dx = 0.25;
fig_y = 0.3;
fig_dy = fig_dx * (1920/1080);
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume perfect value will be plotted
if nargin < 6
    plot_pop = false;
end

% Extract results
s_mat = simResults.s_mat;
h_mat = simResults.h_mat;
l2_mat = log10(simResults.l2_mat);
s_true = simResults.s_true;
h_true = simResults.h_true;
% Prepare plot
Nx = length(x_txts);

% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]); hold on;
% Plot the true value first
plot([-1,Nx+1],[s_true s_true], 'k','LineWidth', 1.5);
% Trim down the data to only the 95% confidence intervals
s_mat = trimToCI(s_mat,keep_frac);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(s_mat,'Whisker',Inf);
% Add colours to the plot
for ii = 1:4
    structfun(@(x) set(x(ii,:), 'color', plot_clrs(ii,:), 'markeredgecolor', plot_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_pop
    s_perfect = simResults.s_perfect;
    hold on;
    for ii = 1:4
        plot(Nx-0.25+0.1*ii,s_perfect(ii),'.','MarkerSize',40, 'MarkerEdgeColor',plot_clrs(ii,:));
    end
end
% Clean up plot
pmin = min(s_mat(:));
pmax = max(s_mat(:));
ymin = max([y_absmin; pmin-0.1*(pmax-pmin)]);
ymax = min([y_absmax; pmax+0.1*(pmax-pmin)]);
axis([0.5 Nx+0.5 ymin ymax]);
xticks(1:Nx);
xticklabels(x_txts);
set(gca,'FontSize',20,'LineWidth',2);
if plot_pop
    xlabel('$N$','FontSize',24,'Interpreter','latex')
end
ylabel('$s$','Fontsize',24,'Interpreter','latex');
title(['Estimated $s$ (',title_txt,')'],'FontSize',20,'Interpreter','latex');

% Save the figure and then close it
if save_figs
    saveas(fig_obj, [fig_filenames{1},'.eps'], 'epsc');
    close(fig_obj);
end

% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]); hold on;
% Plot the true value first
plot([-1,Nx+1],[h_true h_true], 'k','LineWidth', 1.5);
% Trim down the data to only the 95% confidence intervals
h_mat = trimToCI(h_mat,keep_frac);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(h_mat,'Whisker',Inf);
% Add colours to the plot
for ii = 1:4
    structfun(@(x) set(x(ii,:), 'color', plot_clrs(ii,:), 'markeredgecolor', plot_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_pop
    h_perfect = simResults.h_perfect;
    hold on;
    for ii = 1:4
        plot(Nx-0.25+0.1*ii,h_perfect(ii),'.','MarkerSize',40, 'MarkerEdgeColor',plot_clrs(ii,:));
    end
end
plot([-1,Nx+1],[h_true h_true], 'k','LineWidth', 1.5);
% Clean up plot
pmin = min(h_mat(:));
pmax = max(h_mat(:));
ymin = max([y_absmin; pmin-0.1*(pmax-pmin)]);
ymax = min([y_absmax; pmax+0.1*(pmax-pmin)]);
axis([0.5 Nx+0.5 ymin ymax]);
xticks(1:Nx);
xticklabels(x_txts);
set(gca,'FontSize',20,'LineWidth',2);
if plot_pop
    xlabel('$N$','FontSize',24,'Interpreter','latex')
end
ylabel('$h$','Fontsize',24,'Interpreter','latex');
title(['Estimated $h$ (',title_txt,')'],'FontSize',20,'Interpreter','latex');

% Save the figure and then close it
if save_figs
    saveas(fig_obj, [fig_filenames{2},'.eps'], 'epsc');
    close(fig_obj);
end

% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]); hold on;
% Trim down the data to only the 95% confidence intervals
l2_mat = trimToCI(l2_mat,keep_frac);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(l2_mat,'Whisker',Inf);
% Add colours to the plot
for ii = 1:4
    structfun(@(x) set(x(ii,:), 'color', plot_clrs(ii,:), 'markeredgecolor', plot_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_pop
    l2_perfect = log10(simResults.l2_perfect);
    hold on;
    for ii = 1:4
        plot(Nx-0.25+0.1*ii,l2_perfect(ii),'.','MarkerSize',40, 'MarkerEdgeColor',plot_clrs(ii,:));
    end
    axis([0.5 Nx+0.5  -0.1+min([l2_mat(:);l2_perfect(:)]) 0.1+max([l2_mat(:); l2_perfect(:)]) ]);
else
    axis([0.5 Nx+0.5  -0.1+min(l2_mat(:)) 0.1+max(l2_mat(:)) ]);
end

% Clean up plot
xticks(1:Nx);
xticklabels(x_txts);
set(gca,'FontSize',20,'LineWidth',2);
if plot_pop
    xlabel('$N$','FontSize',24,'Interpreter','latex')
end
ylabel('$\log_{10} \mathrm{RMSE}$','Fontsize',24,'Interpreter','latex');
title(['Trajectory $L_2$ Error (',title_txt,')'],'FontSize',20,'Interpreter','latex');

% Save the figure and then close it
if save_figs
    saveas(fig_obj, [fig_filenames{3},'.eps'], 'epsc');
    close(fig_obj);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [s,h] = extract2AlleleParams(A,type)
%
% This function extracts the parameters (s,h) from the payoff matrix that
% was provided as input, using the invariant rules allowed by this type
% of selection are used to normalise and (if possible) symmeterise the 
% matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalise and symmeterise the matrix
if type == 1
    A = A - A(2,2) + 1;
    A(:,1) = A(:,1) - A(2,1) + A(1,2);
else
    A = A / A(2,2);
    A = (A+A')/2;             % Forcibly symmeterise for type II
end
            
% Extract the values of s and h 
s = A(1,1) - 1;
h = (A(1,2) - 1) / s;

