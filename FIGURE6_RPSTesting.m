function FIGURE6_RPSTesting(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying selection strength and
% dominance for a bi-allelic inheritance scenario

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify whether to save figures or just leave them all displayed
save_figs = false;

% Specify the population level to compare methods for
pop_level = 2;

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INPUT HANDLING

% Assume not regenerating data unless specifically requested
if nargin < 1
    regenerate = false;
%elseif regenerate
%    userYN = input("All data will be regenerated. Please confirm Y/N: ",'s');
%    if strcmpi(userYN,'y')
%        regenerate = true;
%    else
%        regenerate = false;
%    end
end


%%% PROBLEM PREPARATION

% Load in the rock-paper-scissors problems and store these in a list of
% problems to run
load('problems_RPS.mat','RPS_balanced','RPS_attract','RPS_repel');
problem_list = {RPS_balanced, RPS_attract, RPS_repel};

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2] };

% Provide the filename prefix
data_prefix = 'DATA_RPS';
% Provide the filenames associated with each problem in the list
problem_filenames = {'BAL','ATR','REP'};
% Provide the additional text added to mark which method was used
method_filenames = {'orders01','orders1','orders012'};
% Title text for the methods and problems
problem_titles = {'Balanced', 'Attracting', 'Repelling'};
method_titles = {'$p \leq 1$', '$p = 1$', '$p \leq 2$'};

% Count provided number of methods
N_methods = length(method_filenames);


%%% GENERATE DATA IF IT IS NOT PRESENT (OR IF ASKED)

% Loop over each problem in the problems list
for p = 1:length(problem_list)
    
    % Read out the current problem from the list
    problem = problem_list{p};
        
    % Loop over each method to run and handle each independently
    for m = 1:N_methods
    
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
    
        % Specify options for this method
        ELoptions.library_orders = library_orders{m};
        ELoptions.symmetric_payoff = false;
    
        % Generate data for this problem if data not present or requested to
        if ~exist(this_filename,'file') || regenerate
            simResults = runSims( problem, ELoptions, @(K,Ftexts,type) extractRPSParams(fitnessToPayoff(K,Ftexts,3), type), pop_sizes, N_rep );
            save(this_filename, 'simResults', 'ELoptions');
        end
        
    end
end

% Safety pause
pause(1);


%%% PLOT FIGURE 3A: TRAJECTORY FIT QUALITY

% Prepare x-ticks for all figures (method descriptions)
x_txts = {'p \leq 1', 'p = 1', 'p \leq 2'};

% Loop over each problem in the problems list
for p = 1:length(problem_list)
        
    % Load in the data for each method
    for m = 1:N_methods
        
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
        % Load data
        load(this_filename, 'simResults');
        
        % Use the first dataset to prepare storage and check population size
        if m == 1
            
            % Count number of replicates and warn user if different to above
            N_reps_here = size(simResults.l2_mat,3);
            if N_rep ~= N_reps_here
                warning('Number of replicates in data does not match number of replicates specified in data generator. Check data and re-generate if required.');
                N_rep = N_reps_here;
            end
            
            % Prepare results storage
            l2_mat = zeros(N_methods, 4, N_rep);
            
            % Note population size
            N_pop = simResults.pop_sizes(pop_level);
            
        end
        
        % Grab out the data for the requested population level and store
        l2_mat(m, :, :) = simResults.l2_mat(pop_level,:,:);
        
    end
    
    % Prepare the title for the figure for this problem
    title_txt = [problem_titles{p},', N = ',num2str(N_pop)];

    % Plot the error boxplots for this problem
    plotRMSEBoxplot( struct('l2_mat', l2_mat), x_txts, title_txt, ['FIG6a_',problem_filenames{p},'_methodRMSE'], save_figs, false);
    
end



%%% PLOT FIGURE 3B: RECOVERY OF PARAMETERS FOR DIFFERENT POPULATION SIZES

% Prepare the x tick labels
for n = 1:length(simResults.pop_sizes)
    x_txts{n} = num2str(simResults.pop_sizes(n));
end
x_txts{end+1} = '\infty';

% Loop across all problems
for p = 1:length(problem_list)

    % Extract the true parameters for this problem
    A = problem.fitness(eye(3));
    params_true = extractRPSParams(A,1);
    
    % Loop across all methods figures
    for m = 1:N_methods

        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
        % Load data
        load(this_filename, 'simResults');
        % Add the true parameter values to the results struct
        simResults.params_true = params_true;
        % Carry out the plot
        plotParamsBoxplot( simResults, x_txts, [problem_titles{p},', ',method_titles{m}], ['FIG6b_',problem_filenames{p},'_pop'], save_figs, true );
        plotRMSEBoxplot( simResults, x_txts, title_txt, ['FIG6c_',problem_filenames{p},'_popRMSE'], save_figs, false);
        
    end 
end