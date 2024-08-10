function FIGURE_basicAlleleTesting(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying selection strength and
% dominance for a bi-allelic inheritance scenario

% Specify a random seed for consistency
rng(7);

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 1000;

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INPUT HANDLING

% Assume not regenerating data unless specifically requested
if nargin < 1
    regenerate = false;
elseif regenerate
    userYN = input("All data will be regenerated. Please confirm Y/N: ");
    if ~STRCMP(userYN,'y')
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


%%% RUNNING PROBLEMS UNDER DIFFERENT SETTINGS

% General linear equation learning
title_txt = 'Up to First Order';
filename = 'DATA_Npop_orders01';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.orders = [0,1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotOverPopSizes(simResults, title_txt);

% Strictly payoff matrix elements
title_txt = 'Strictly First Order';
filename = 'DATA_Npop_orders1';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.orders = [1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotOverPopSizes(simResults, title_txt);

% Quadratic fitness function
title_txt = 'Up to Second Order';
filename = 'DATA_Npop_orders012';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = false;
    ELoptions.orders = [0,1,2];
    simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotOverPopSizes(simResults, title_txt);

% Specifically a symmetric payoff matrix
title_txt = 'Symmetric Payoff';
filename = 'DATA_Npop_symm';
% Regenerate if requested or data not present
if ~exist([filename,'.mat'],'file') || regenerate
    ELoptions.symmetric_payoff = true;
    ELoptions.orders = [1];
    simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
else
    load([filename,'.mat'],'simResults');
end    
plotOverPopSizes(simResults, title_txt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename)
%
% This subfunction generates data for equation learning applied to the
% provided allele inheritance problems, returning the estimated values for
% the parameters (s,h) for different population sizes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over the list of population sizes, learning parameters s and h for
% each
N_pops = length(pop_sizes);
s_mat = zeros(N_pops, 2, N_rep);
h_mat = zeros(N_pops, 2, N_rep);
for k = 1:N_pops
   
    % Loop over the two types of selection and generate results for each
    for type = [1,2]
    
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
        
            % Extract learned payoff matrix and normalise using a_22 = 1
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
             
        end 
    end    
    
end

% Now generate "perfect" replicator data with known derivative and estimate
% using it for reference
s_perfect = zeros(1,2);
h_perfect = zeros(1,2);
for type = [1,2]
    
    % Load in problem
    problem = problems{type};
    
    % Generate the replicator data, including derivatives
    traj = generateTrajectories(problem);
    X_true = traj.rep.X;
    Xdash_true = zeros(size(X_true));
    for j = 1:size(X_true,2)
        Xdash_true(:,j) = replicatorRHS(X_true(:,j),problem.fitness,type);
    end
    
    ELoptions.Xdash_data = Xdash_true;
    [fitfun, Flib, K, texts] = evolutionLearning(traj.rep, type, ELoptions);
    
    % Extract learned payoff matrix and normalise using a_22 = 1
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
simResults = struct('s_mat',s_mat,'h_mat',h_mat,'s_perfect',s_perfect,'h_perfect',h_perfect,'s_true',s_true,'h_true',h_true,'pop_sizes',pop_sizes);
save([filename,'.mat'],'simResults');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function plotOverPopSizes(simResults, title_txt)
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
ymin = -0.75;
ymax = 1.75;
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract results
s_mat = simResults.s_mat;
h_mat = simResults.h_mat;
s_perfect = simResults.s_perfect;
h_perfect = simResults.h_perfect;
s_true = simResults.s_true;
h_true = simResults.h_true;
pop_sizes = simResults.pop_sizes;
% Prepare plot
N_pops = length(pop_sizes);

% Prepare the labels for each population size
x_txts = cell(1,N_pops);
for k = 1:N_pops
    x_txts{k} = ['N = ',num2str(pop_sizes(k))];
end
x_txts{N_pops+1} = 'Perfect Data';

% Prepare the figure
figure; hold on;
% Plot the true value first
plot([-1,N_pops+2],[s_true s_true], 'k','LineWidth', 1.5);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(s_mat);
% Add colours to the plot
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:)), boxplot_obj);
end
% Add to the plot the parameter learned for perfect data
hold on;
plot(N_pops+0.9,s_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
plot(N_pops+1.1,s_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
% Clean up plot
axis([0 N_pops+1.2 ymin ymax]);
xticks(1:N_pops+1);
xticklabels(x_txts);
set(gca,'FontSize',20);
ylabel('s','Fontsize',24);
title(['Estimation of Selective Advantage --- ',title_txt],'FontSize',20);

% Prepare the figure
figure; hold on;
% Plot the true value first
plot([-1,N_pops+2],[h_true h_true], 'k','LineWidth', 1.5);
% Plot the boxplot showing main data
boxplot_obj = boxplot2(h_mat);
% Add colours to the plot
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:)), boxplot_obj);
end
hold on;
plot(N_pops+0.9,h_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
plot(N_pops+1.1,h_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
plot([-1,N_pops+2],[h_true h_true], 'k','LineWidth', 1.5);
% Clean up plot
axis([0 N_pops+1.2 ymin ymax]);
xticks(1:N_pops+1);
xticklabels(x_txts);
set(gca,'FontSize',20);
ylabel('h','Fontsize',24);
title(['Estimation of Dominance --- ',title_txt],'FontSize',20);