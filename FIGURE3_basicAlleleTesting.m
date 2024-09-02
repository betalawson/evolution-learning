function FIGURE3_basicAlleleTesting(regenerate)
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
load('problems_basic.mat','allele_basic1','allele_basic2');

% Store type I and type II in a single cell array
problems = {allele_basic1, allele_basic2};


%%% EFFECTS OF POPULATION SIZE (ACROSS METHODS)

% General linear equation learning
title_txt = 'Up to First Order';
filename = 'DATA_Npop_orders01';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.library_orders = [0,1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end

% This is the first problem so use it to set up the population counts
N_pops = length(simResults.pop_sizes);
x_txts = cell(1,N_pops);
for k = 1:N_pops
    x_txts{k} = ['N = 10^{',num2str(log10(simResults.pop_sizes(k))),'}'];
end
x_txts{N_pops+1} = 'Perfect Data';

% Actually plot
plotBoxPlot(simResults, x_txts, title_txt);

% Strictly payoff matrix elements
title_txt = 'Strictly First Order';
filename = 'DATA_Npop_orders1';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.library_orders = [1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotBoxPlot(simResults, x_txts, title_txt);

% Quadratic fitness function
title_txt = 'Up to Second Order';
filename = 'DATA_Npop_orders012';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.library_orders = [0,1,2];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotBoxPlot(simResults, x_txts, title_txt);

% Specifically a symmetric payoff matrix
title_txt = 'Symmetric Payoff';
filename = 'DATA_Npop_symm';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = true;
    ELoptions.library_orders = [1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, test_types, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotBoxPlot(simResults, x_txts, title_txt);


%%% EFFECTS OF METHOD CHOICE (FOR A SINGLE POPULATION SIZE)

% Specify the population level to compare methods for
pop_level = 2;

% Define the list of filenames to extract data from
filenames = {'DATA_Npop_orders1', 'DATA_Npop_orders01', 'DATA_Npop_orders012', 'DATA_Npop_symm'};

% Define the methods text used on the figure
method_txts = {'1st Order', '0th-1st Order', '0th-2nd order', 'Symmetric'};

% Load in the first dataset to set up the structure
load([filenames{1},'.mat'],'simResults');
s_mat = zeros([length(filenames), 2, size(simResults.s_mat,3)]);
h_mat = zeros([length(filenames), 2, size(simResults.s_mat,3)]);
l2_mat = zeros([length(filenames), 2, size(simResults.s_mat,3)]);
s_mat(1,:,:) = simResults.s_mat(pop_level,:,:);
h_mat(1,:,:) = simResults.h_mat(pop_level,:,:);
l2_mat(1,:,:) = simResults.l2_mat(pop_level,:,:);

% Load in each other dataset and extract results for this population level
for m = 2:length(filenames)
    load([filenames{m},'.mat'],'simResults');
    s_mat(m, :, :) = simResults.s_mat(pop_level,:,:);
    h_mat(m, :, :) = simResults.h_mat(pop_level,:,:);
    l2_mat(m, :, :) = simResults.l2_mat(pop_level,:,:);
end

% Plot these results
plotBoxPlot(struct('s_mat',s_mat,'h_mat',h_mat,'l2_mat',l2_mat,'s_true',simResults.s_true,'h_true',simResults.h_true), method_txts, ['N = ',num2str(simResults.pop_sizes(pop_level))], false);



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
s_mat = zeros(N_pops, 2, N_rep);
h_mat = zeros(N_pops, 2, N_rep);
l2_mat = zeros(N_pops, 2, N_rep);
for k = 1:N_pops
   
    % Loop over the two types of selection and generate results for each
    for type = test_types
    
        % Overwrite the population size for this problem
        problem = problems{type};
        problem.N_pop = pop_sizes(k);
    
        % Repeatedly generate new Wright-Fisher data and learn s and h
        parfor r = 1:N_rep
        
            % Generate the Wright-Fisher data
            traj = generateTrajectories(problem);
        
            % Convert to non-dimensionalised time for equation learning
            WF_nonD = struct('t',traj.WF.t / problem.t_gen, 'X', traj.WF.X);

            % Perform equation learning
            [fitfun, Flib, K, texts] = evolutionLearning(WF_nonD, type, ELoptions);
            
            % Extract payoff matrix and normalise using invariant operations
            payoff = fitnessToPayoff(K,texts,2);            
            if type == 1
                payoff = payoff - payoff(2,2) + 1;
                payoff(:,1) = payoff(:,1) - payoff(2,1) + payoff(1,2);
            else
                payoff = payoff / payoff(2,2);
            end
            
            % Extract the values of s and h and store them 
            s_mat(k,type,r) = payoff(1,1) - 1;
            h_mat(k,type,r) = (0.5 * ( payoff(1,2) + payoff(2,1) ) - 1) / s_mat(k,type,r);
            
            % Simulate the learned dynamics and evaluate the L2 error
            % between the learned dynamics and the actual replicator
            learned_problem = problem;
            learned_problem.fitness = @(X) payoff*X;
            learned_traj = generateTrajectories(learned_problem);
            l2_mat(k,type,r) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));                                    
             
        end 
    end    
    
end

% Now generate "perfect" replicator data with known derivative and estimate
% using it for reference
s_perfect = zeros(1,2);
h_perfect = zeros(1,2);
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
    
    % Provide the perfect data and use it for evolution learning
    ELoptions.Xdash_data = Xdash_true;
    [fitfun, Flib, K, texts] = evolutionLearning(traj.rep, type, ELoptions);
    
    % Extract payoff matrix and normalise using invariant operations
    payoff = fitnessToPayoff(K,texts,2);
    if type == 1
        payoff = payoff - payoff(2,2) + 1;
        payoff(:,1) = payoff(:,1) - payoff(2,1) + payoff(1,2);
    else
        payoff = payoff / payoff(2,2);
    end
    % Extract the values of s and h and store them 
    s_perfect(type) = payoff(1,1) - 1;
    h_perfect(type) = (0.5 * ( payoff(1,2) + payoff(2,1) ) - 1) / s_perfect(type);

    % Find the L2 error to the true selective trajectory
    learned_problem = problem;
    learned_problem.fitness = @(X) payoff*X;
    learned_traj = generateTrajectories(learned_problem);
    l2_perfect(type) = rms(traj.rep.X(1,:) - learned_traj.rep.X(1,:));
    
end

% Extract the true values of s and h from the problem
payoff = problems{1}.fitness(eye(length(problem.X0)));    
if type == 1
    payoff = payoff - payoff(2,2) + 1;
    payoff(:,1) = payoff(:,1) - payoff(2,1) + payoff(1,2);
else
    payoff = payoff / payoff(2,2);
end
% Extract the values of s and h and store them 
s_true = payoff(1,1) - 1;
h_true = (0.5 * ( payoff(1,2) + payoff(2,1) ) - 1) / s_true;

% Store all data in a single results object and save it
simResults = struct('s_mat',s_mat,'h_mat',h_mat,'l2_mat',l2_mat,'s_perfect',s_perfect,'h_perfect',h_perfect,'l2_perfect',l2_perfect,'s_true',s_true,'h_true',h_true,'pop_sizes',pop_sizes);
save([filename,'.mat'],'simResults');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function plotBoxPlot(simResults, x_txts, title_txt, plot_perfect)
%
% This subfunction just loops over the provided population sizes and
% calculates the equation-learned values of the allele inheritance
% parameters across the specified number of replicate realisations of
% Wright-Fisher
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the colours to use in plotting the two selection types
type_clrs = [ [ 1.0, 0.4, 0.4 ];
              [ 0.4, 0.4, 1.0 ]  ];
          
% Specify the minimum and maximum allowable values for s and h on the plots
pmin = -0.75;
pmax = 1.75;
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume perfect value will be plotted
if nargin < 4
    plot_perfect = true;
end

% Extract results
s_mat = simResults.s_mat;
h_mat = simResults.h_mat;
l2_mat = log(simResults.l2_mat);
s_true = simResults.s_true;
h_true = simResults.h_true;
% Prepare plot
Nx = length(x_txts);

% Prepare the figure
figure; hold on;
% Plot the true value first
plot([-1,Nx+1],[s_true s_true], 'k','LineWidth', 1.5);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(s_mat);
% Add colours to the plot
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_perfect
    s_perfect = simResults.s_perfect;
    hold on;
    plot(Nx-0.1,s_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
    plot(Nx+0.1,s_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
end
% Clean up plot
smin = min(s_mat(:));
smax = max(s_mat(:));
ymin = max([pmin; smin-0.1*(smax-smin)]);
ymax = min([pmax; smax+0.1*(smax-smin)]);
axis([0.5 Nx+0.5 ymin ymax]);
xticks(1:Nx);
xticklabels(x_txts);
set(gca,'FontSize',20);
ylabel('s','Fontsize',24);
title(['Selective Advantage (',title_txt,')'],'FontSize',20);

% Prepare the figure
figure; hold on;
% Plot the true value first
plot([-1,Nx+1],[h_true h_true], 'k','LineWidth', 1.5);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(h_mat);
% Add colours to the plot
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_perfect
    h_perfect = simResults.h_perfect;
    hold on;
    plot(Nx-0.1,h_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
    plot(Nx+0.1,h_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
end
plot([-1,Nx+1],[h_true h_true], 'k','LineWidth', 1.5);
% Clean up plot
hmin = min(h_mat(:));
hmax = max(h_mat(:));
ymin = max([pmin; hmin-0.1*(hmax-hmin)]);
ymax = min([pmax; hmax+0.1*(hmax-hmin)]);
axis([0.5 Nx+0.5 ymin ymax]);
xticks(1:Nx);
xticklabels(x_txts);
set(gca,'FontSize',20);
ylabel('h','Fontsize',24);
title(['Dominance (',title_txt,')'],'FontSize',20);

% Prepare the figure
figure; hold on;
% Plot the boxplot showing main data
boxplot_obj = boxplot2(l2_mat);
% Add colours to the plot
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data (if requested)
if plot_perfect
    l2_perfect = log(simResults.l2_perfect);
    hold on;
    plot(Nx-0.1,l2_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
    plot(Nx+0.1,l2_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
    axis([0.5 Nx+0.5  -0.1+min([l2_mat(:);l2_perfect(:)]) 0.1+max([l2_mat(:); l2_perfect(:)]) ]);
else
    axis([0.5 Nx+0.5  -0.1+min(l2_mat(:)) 0.1+max(l2_mat(:)) ]);
end

% Clean up plot
xticks(1:Nx);
xticklabels(x_txts);
set(gca,'FontSize',20);
ylabel('L_2','Fontsize',24);
title(['Trajectory L_2 Error (',title_txt,')'],'FontSize',20);