function FIGURE4_triAlleleTesting(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying selection strength and
% dominance for a bi-allelic inheritance scenario

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify which types of selection to test (1 or 2 or both)
test_types = [1,2];

% Specify whether to save figures or just leave them all displayed
save_figs = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INPUT HANDLING

% Assume not regenerating data unless specifically requested
if nargin < 1
    regenerate = false;
%elseif regenerate
    %userYN = input("All data will be regenerated. Please confirm Y/N: ",'s');
    %if strcmpi(userYN,'y')
    %    regenerate = true;
    %else
    %    regenerate = false;
    %end
end


%%% PROBLEM PREPARATION

% Load in the basic problem - type I and type II selection versions.
% Settings for these will be adjusted as needed
load('problems_inheritance.mat', 'allele_standard', 'allele_transient', 'allele_persistent');

% Create a cell array of cell arrays that contains the problem list
problems_list = { allele_standard, allele_transient, allele_persistent };

% Specify the filename "suffixes" for each problem in the master list
problem_filenames = {'std', 'tran', 'pers'};
problem_names = {'Standard', 'Transient', 'Persistent'};

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;

% Loop over and run each problem from the problem list separately
for P = 1:length(problems_list)
    
    % Load in current problem from master list
    base_problem = problems_list{P};
    
    % Create a type I and type II selection version of the problem
    problems = cell(1,2);
    for type = 1:2
        problem = base_problem;
        problem.selection_type = type;
        problems{type} = problem;
    end
        
%     % Using method: Up to First Order
%     data_filename = ['DATA_IN',problem_filenames{P},'_Npop_01'];
%     % Regenerate if requested or data not present
%     if ~exist([data_filename,'.mat'],'file') || regenerate
%         ELoptions.symmetric_payoff = false;
%         ELoptions.library_orders = [0,1];
%         simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
%     else
%         load([data_filename,'.mat'],'simResults');
%     end
    
%     % Using method: Up to First Order
%     data_filename = ['DATA_IN',problem_filenames{P},'_Npop_1'];
%     % Regenerate if requested or data not present
%     if ~exist([data_filename,'.mat'],'file') || regenerate
%         ELoptions.symmetric_payoff = false;
%         ELoptions.library_orders = [1];
%         simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
%     else
%         load([data_filename,'.mat'],'simResults');
%     end   
    
%     % Using method: Up to First Order
%     data_filename = ['DATA_IN',problem_filenames{P},'_Npop_012'];
%     % Regenerate if requested or data not present
%     if ~exist([data_filename,'.mat'],'file') || regenerate
%         ELoptions.symmetric_payoff = false;
%         ELoptions.library_orders = [0,1,2];
%         simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
%     else
%         load([data_filename,'.mat'],'simResults');
%     end   
    
    % Using method: Symmetric payoff matrix
    data_filename = ['DATA_IN',problem_filenames{P},'_Npop_symm'];
    % Regenerate if requested or data not present
    if ~exist([data_filename,'.mat'],'file') || regenerate
        ELoptions.symmetric_payoff = true;
        ELoptions.library_orders = [1];
        simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, data_filename);
    else
        load([data_filename,'.mat'],'simResults');
    end    
    
    % Prepare the x labels
    N_pops = length(simResults.pop_sizes);
    x_txts = cell(1,N_pops);
    for k = 1:N_pops
        x_txts{k} = ['10^{',num2str(log10(simResults.pop_sizes(k))),'}'];
    end
    x_txts{N_pops+1} = '\infty';
    
    % List the figure filenames
    param_filenames = {'s1','s2','h12','h13','h23','l2'};
    num_suffixes = {'a','b','c','d','e','f'};
    N_figs = length(param_filenames);
    fig_filenames = cell(1,N_figs);
    for k = 1:N_figs
        fig_filenames{k} = ['FIG4',num_suffixes{k},'_triallele_',param_filenames{k}];
    end
    title_txt = problem_names{P};
    plotBoxPlot(simResults, x_txts, title_txt, fig_filenames, save_figs, true)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, filename)
%
% This subfunction generates data for equation learning applied to the
% provided allele inheritance problems, returning the estimated values for
% the parameters s_1, s_2, h_12, h_13, h_23 for the different population
% sizes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify a random seed for consistency - this should give the same
% trajectories for each different approach to match, also
rng(7);

% Loop over the list of population sizes, learning parameters for each
N_pops = length(pop_sizes);
s1LS_mat = zeros(N_pops, 2, N_rep);
s2LS_mat = zeros(N_pops, 2, N_rep);
h12LS_mat = zeros(N_pops, 2, N_rep);
h13LS_mat = zeros(N_pops, 2, N_rep);
h23LS_mat = zeros(N_pops, 2, N_rep);
l2LS_mat = zeros(N_pops, 2, N_rep);
s1EL_mat = zeros(N_pops, 2, N_rep);
s2EL_mat = zeros(N_pops, 2, N_rep);
h12EL_mat = zeros(N_pops, 2, N_rep);
h13EL_mat = zeros(N_pops, 2, N_rep);
h23EL_mat = zeros(N_pops, 2, N_rep);
l2EL_mat = zeros(N_pops, 2, N_rep);
LS_time_tot = 0;
EL_time_tot = 0;
for k = 1:N_pops
   
    % Loop over the two types of selection and generate results for each
    for type = test_types
    
        % Overwrite the population size for this problem
        problem = problems{type};
        problem.N_pop = pop_sizes(k);
    
        % Initialise time storage
        LS_times = zeros(1,N_rep);
        EL_times = zeros(1,N_rep);
        
        % Repeatedly generate new Wright-Fisher data and learn s and h
        parfor r = 1:N_rep
        
            % Generate the Wright-Fisher and replicator data
            traj = generateTrajectories(problem);
        
            % Convert to non-dimensionalised time for equation learning
            WF_nonD = struct('t',traj.WF.t / problem.t_gen, 'X', traj.WF.X);

            % Perform equation learning
            tic;
            [Flib, K, texts] = LSevolutionLearning(WF_nonD, type, ELoptions);
            LS_times(r) = toc;
            
            % Extract payoff matrix and normalise using invariant operations
            A = fitnessToPayoff(K,texts,3);
            [s1,s2,h12,h13,h23] = extract3AlleleParams(A,type);

            % Extract the values of s and h and store them 
            s1LS_mat(k,type,r) = s1;
            s2LS_mat(k,type,r) = s2;
            h12LS_mat(k,type,r) = h12;
            h13LS_mat(k,type,r) = h13;
            h23LS_mat(k,type,r) = h23;
            
            % Simulate the learned dynamics and evaluate the L2 error
            % between the learned dynamics and the actual replicator
            learned_problem = problem;
            learned_problem.fitness = @(X) libraryRHS(X,Flib,K);
            learned_traj = generateTrajectories(learned_problem,'rep');
            l2LS_mat(k,type,r) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));   
            
            % Now use our gradient matching method
            tic;
            [fitfun, Flib, K, texts] = evolutionLearning(WF_nonD, type, ELoptions);
            EL_times(r) = toc;
            
            % Extract payoff matrix and normalise using invariant operations
            A = fitnessToPayoff(K,texts,3);
            [s1,s2,h12,h13,h23] = extract3AlleleParams(A,type);

            % Extract the values of s and h and store them 
            s1EL_mat(k,type,r) = s1;
            s2EL_mat(k,type,r) = s2;
            h12EL_mat(k,type,r) = h12;
            h13EL_mat(k,type,r) = h13;
            h23EL_mat(k,type,r) = h23;
            
            % Simulate the learned dynamics and evaluate the L2 error
            % between the learned dynamics and the actual replicator
            learned_problem = problem;
            learned_problem.fitness = @(X) libraryRHS(X,Flib,K);
            learned_traj = generateTrajectories(learned_problem,'rep');
            l2EL_mat(k,type,r) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));                       
             
        end 
        
        LS_time_tot = LS_time_tot + sum(LS_times);
        EL_time_tot = EL_time_tot + sum(EL_times);
        
    end    
    
end

% Now gather together the data
s1_mat(:,1:2,:) = s1LS_mat;
s1_mat(:,3:4,:) = s1EL_mat;
s2_mat(:,1:2,:) = s2LS_mat;
s2_mat(:,3:4,:) = s2EL_mat;
h12_mat(:,1:2,:) = h12LS_mat;
h12_mat(:,3:4,:) = h12EL_mat;
h13_mat(:,1:2,:) = h13LS_mat;
h13_mat(:,3:4,:) = h13EL_mat;
h23_mat(:,1:2,:) = h23LS_mat;
h23_mat(:,3:4,:) = h23EL_mat;
l2_mat(:,1:2,:) = l2LS_mat;
l2_mat(:,3:4,:) = l2EL_mat;

% Now generate "perfect" replicator data with known derivative and estimate
% using it for reference
s1_perfect = zeros(1,2);
s2_perfect = zeros(1,2);
h12_perfect = zeros(1,2);
h13_perfect = zeros(1,2);
h23_perfect = zeros(1,2);
l2_perfect = zeros(1,2);
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
    A = fitnessToPayoff(K,texts,3);
    [s1,s2,h12,h13,h23] = extract3AlleleParams(A,type);
    
    % Extract the values of s and h and store them 
    s1_perfect(type) = s1;
    s2_perfect(type) = s2;
    h12_perfect(type) = h12;
    h13_perfect(type) = h13;
    h23_perfect(type) = h23;
    
    % Find the L2 error to the true selective trajectory
    learned_problem = problem;
    learned_problem.fitness = @(X) libraryRHS(X,Flib, K);
    learned_traj = generateTrajectories(learned_problem);
    l2_perfect(type) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));
    
    % Provide the perfect data and use it for evolution learning
    ELoptions.Xdash_data = Xdash_true;
    [fitfun, Flib, K, texts] = evolutionLearning(traj.rep, type, ELoptions);
    
    % Extract payoff matrix and normalise using invariant operations
    A = fitnessToPayoff(K,texts,3);
    [s1,s2,h12,h13,h23] = extract3AlleleParams(A,type);
    
    % Extract the values of s and h and store them 
    s1_perfect(type+2) = s1;
    s2_perfect(type+2) = s2;
    h12_perfect(type+2) = h12;
    h13_perfect(type+2) = h13;
    h23_perfect(type+2) = h23;

    % Find the L2 error to the true selective trajectory
    learned_problem = problem;
    learned_problem.fitness = @(X) libraryRHS(X,Flib, K);
    learned_traj = generateTrajectories(learned_problem);
    l2_perfect(type+2) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));
    
end

% Extract the true values of s and h from the problem
A = problem.fitness(eye(length(problem.X0)));
[s1_true,s2_true,h12_true,h13_true,h23_true] = extract3AlleleParams(A,type);

% Store all data in a single results object and save it
simResults = struct('s1_mat',s1_mat,'s2_mat',s2_mat,'h12_mat',h12_mat,'h13_mat',h13_mat,'h23_mat',h23_mat,'l2_mat',l2_mat,'s1_perfect',s1_perfect,'s2_perfect',s2_perfect,'h12_perfect',h12_perfect,'h13_perfect',h13_perfect,'h23_perfect',h23_perfect,'l2_perfect',l2_perfect,'s1_true',s1_true,'s2_true',s2_true,'h12_true',h12_true,'h13_true',h13_true,'h23_true',h23_true,'pop_sizes',pop_sizes,'LS_time_tot',LS_time_tot,'EL_time_tot',EL_time_tot);
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
type_clrs = [   [ 0.6, 0.0, 0.0 ];          % Dark red
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
s1_mat = trimToCI(simResults.s1_mat,keep_frac);
s2_mat = trimToCI(simResults.s2_mat,keep_frac);
h12_mat = trimToCI(simResults.h12_mat,keep_frac);
h13_mat = trimToCI(simResults.h13_mat,keep_frac);
h23_mat = trimToCI(simResults.h23_mat,keep_frac);
l2_mat = trimToCI(log10(simResults.l2_mat),keep_frac);
s1_perfect = simResults.s1_perfect;
s2_perfect = simResults.s2_perfect;
h12_perfect = simResults.h12_perfect;
h13_perfect = simResults.h13_perfect;
h23_perfect = simResults.h23_perfect;
l2_perfect = log10(simResults.l2_perfect);
s1_true = simResults.s1_true;
s2_true = simResults.s2_true;
h12_true = simResults.h12_true;
h13_true = simResults.h13_true;
h23_true = simResults.h23_true;

% Gather data together into cell arrays that enable the below for loop
p_mats = {s1_mat, s2_mat, h12_mat, h13_mat, h23_mat};
p_trues = {s1_true, s2_true, h12_true, h13_true, h23_true};
p_perfects = {s1_perfect, s2_perfect, h12_perfect, h13_perfect, h23_perfect};

% Specify the parameter names (letters and descriptions)
%%% (DESCRIPTIONS CURRNETLY UNUSED, TOO LONG FOR PLOTS )
p_names = {'$s_1$', '$s_2$', '$h_{12}$', '$h_{13}$', '$h_{23}$'};
p_descs = {'Selective Advantage of Allele 1', 'Selective Advantage of Allele 2', 'Dominance of Allele 1 over 2', 'Dominance of Allele 1 over 3', 'Dominance of Allele 2 over 3'}; 

% Prepare plot
Nx = length(x_txts);

% Loop over each parameter in the cell arrays
for pnum = 1:length(p_mats)

    % Extract this parameter's data
    p_mat = p_mats{pnum};
    p_true = p_trues{pnum};
    p_perfect = p_perfects{pnum};
    
    % Prepare the figure
    fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]); hold on;
    % Plot the true value first
    plot([-1,Nx+1],[p_true p_true], 'k','LineWidth', 1.5);
    % Plot the boxplot showing main data
    boxplot_obj = boxplot2(p_mat,'Whisker',Inf);
        
    % Add colours to the plot
    for ii = 1:4
        structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
    end
    
    % Add to the plot the parameter learned for perfect data (if requested)
    if plot_pop
        plot(Nx-0.1,p_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
        plot(Nx+0.1,p_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
    end
    
    % Determine appropriate y-axis limits for this parameter
    pmin = min(p_mat(:));
    pmax = max(p_mat(:));
    ymin = max([y_absmin; pmin-0.1*(pmax-pmin)]);
    ymax = min([y_absmax; pmax+0.1*(pmax-pmin)]);
    
    % Clean up plot
    axis([0.5 Nx+0.5 ymin ymax]);
    xticks(1:Nx);
    xticklabels(x_txts);
    set(gca,'FontSize',20,'LineWidth',2);
    if plot_pop
        xlabel('$N$','FontSize',24,'Interpreter','latex')
    end
    ylabel(p_names{pnum},'Fontsize',24,'Interpreter','latex');
    title(['Estimated ',p_names{pnum},' (',title_txt,')'],'FontSize',20,'Interpreter','latex');

    % Save the figure and then close it
    if save_figs
        saveas(fig_obj, [fig_filenames{pnum},'.eps'], 'epsc');
        close(fig_obj);
    end
    
end

% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]); hold on;
% Trim down the data to only the 95% confidence intervals
l2_mat = trimToCI(l2_mat,keep_frac);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(l2_mat,'Whisker',Inf);
% Add colours to the plot
for ii = 1:4
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_pop
    hold on;
    plot(Nx-0.15,l2_perfect(1),'.','MarkerSize',30, 'MarkerEdgeColor',type_clrs(1,:));
    plot(Nx-0.05,l2_perfect(2),'.','MarkerSize',30, 'MarkerEdgeColor',type_clrs(2,:));
    plot(Nx+0.05,l2_perfect(3),'.','MarkerSize',30, 'MarkerEdgeColor',type_clrs(3,:));
    plot(Nx+0.15,l2_perfect(4),'.','MarkerSize',30, 'MarkerEdgeColor',type_clrs(4,:));
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
function [s1,s2,h12,h13,h23] = extract3AlleleParams(A,type)
%
% This function extracts the parameters (s,h) from a provided fitness
% library of polynomial terms. The payoff matrix is first extracted, and
% then the invariant rules allowed by this type of selection are used to
% normalise and (if possible) symmeterise the matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract learned payoff matrix and normalise using a_33 = 1
if type == 1
    
    %%% FIND CLOSEST-TO-SYMMETRIC "SIMILAR" MATRIX
    %
    % A similar matrix here means a matrix that can be reached
    % using shifts to columns, as the dynamics are invariant
    % to such shifts.
    %
    % Symmetry conditions are Cz = b
    %    C - matrix showing which shifts apply where
    %    z - shifts to apply
    %    b - symmetry checks
    
    % Closed form for C
    C = [1 -1 0; 1 0 -1; 0 1 -1];
    % Symmetry checks
    b = [ A(1,2) - A(2,1) ; A(1,3) - A(3,1) ; A(2,3) - A(3,2) ];
    % Find minimum-norm x that gets as close as possible to
    % satisfying symmetry checks
    z = lsqminnorm(C, b);
    
    % Adjust the payoff matrix using these shifts
    A = A + z';
    
    % Now adjust the payoff matrix so the bottom-right element
    % is unity to match the base matrix definition
    A = A - A(3,3) + 1;
    
else
    
    % Type II selection only allows constant scaling, so we
    % cannot do anything to push towards symmetry, only
    % normalise using bottom-right element
    A = A / A(3,3);
    
end

% Symmeterise the matrix
A = (A + A') / 2;

% Extract the values of s and h and store them
s1 = A(1,1) - 1;
s2 = A(2,2) - 1;
h12 = ( A(1,2) - 1 - s2 ) / (s1 - s2);
h13 = ( A(1,3) - 1 ) / s1;
h23 = ( A(2,3) - 1 ) / s2;