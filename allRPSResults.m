function allRPSResults(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying rock paper scissors style
% dynamics, which naturally oscillate

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify whether to save figures or just leave them all displayed
save_figs = false;

% Specify the population level to compare methods for
pop_level = 2;

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

% Specify how many generations to simulate for prediction error checking
sim_gen = 200;

% Dataset name
data_prefix = 'DATA_RPS_';

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

% Load in the rock-paper-scissors problems
load('problems_RPS.mat','RPS_balanced','RPS_attract','RPS_repel');

% Create a cell array of cell arrays that contains the problem list
problems_list = {RPS_balanced, RPS_attract, RPS_repel};

% Specify the filename "suffixes" for each problem in the master list
problem_filenames = {'bal','atr','rep'};
problem_names = {'Balanced', 'Attracting', 'Repelling'};

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2] };

% Provide the additional text added to mark which method was used
method_filenames = {'orders01','orders1','orders012'};

% Count provided number of methods and problems
N_methods = length(method_filenames);
N_problems = length(problems_list);


%%% GENERATE DATA IF IT IS NOT PRESENT (OR IF ASKED)

% Loop over each problem in the problem list
for p = 1:N_problems

    % Loop over each method to run and handle each independently
    for m = 1:N_methods
    
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
    
        % Specify options for this method
        ELoptions.library_orders = library_orders{m};
        ELoptions.symmetric_payoff = false;
        ELoptions.force_positive = true;
    
        % Generate data for this problem if data not present or requested to
        if ~exist(this_filename,'file') || regenerate
            simResults = runSims( problems_list{p}, ELoptions, @(library,type) extractRPSParams(libraryToPayoff(library,3), type), pop_sizes, N_rep, sim_gen );
            save(this_filename, 'simResults', 'ELoptions');
        end
   
    end
end

% Safety pause
pause(1);


%%% PLOT FIGURE 6A: TRAJECTORY FIT QUALITY

% Loop over each problem
for p = 1:N_problems

    % Load in the data for each method
    for m = 1:N_methods
    
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
        % Load data
        load(this_filename, 'simResults');
       
        % Use the first dataset to prepare storage and check population size
        if m == 1
        
            % Count number of replicates and warn user if different to above
            N_reps_here = size(simResults.obsl2_mat,3);
            if N_rep ~= N_reps_here
                warning('Number of replicates in data does not match number of replicates specified in data generator. Check data and re-generate if required.');
                N_rep = N_reps_here;
            end
        
            % Prepare results storage
            obsl2_mat = zeros(N_methods, 4, N_rep);
            predl2_mat = zeros(N_methods, 4, N_rep);
        
            % Note population size
            N_pop = simResults.pop_sizes(pop_level);
        
        end
        
        % Grab out the data for the requested population level and store
        obsl2_mat(m, :, :) = simResults.obsl2_mat(pop_level,:,:);
        predl2_mat(m, :, :) = simResults.predl2_mat(pop_level,:,:);
    
    end
    
    % Prepare all text for figure
    x_txts = {'p \leq 1', 'p = 1', 'p \leq 2'};
    title_txt = [problem_names{p},' N = ',num2str(N_pop)];
    % Plot the error boxplots
    plotRMSEBoxplot( struct('l2_mat', obsl2_mat), x_txts, [title_txt, ' (Observed)'], ['FIG6a_',problem_filenames{p},'_observedRMSE'], save_figs, false);
    plotRMSEBoxplot( struct('l2_mat', predl2_mat), x_txts, [title_txt, ' (Forecasted)'], ['FIG6a_',problem_filenames{p},'_forecastRMSE'], save_figs, false);
    
end


%%% PLOTS FOR PARAMETER RECOVERY AND ERROR

% Provide text for methods
method_txts = {'$p \leq 1$', '$p = 1$', '$p \leq 2$'};

% Loop over all problems
for p = 1:N_problems
    
    % Get the true values of the parameters
    A = problems_list{p}.fitness(eye(3));
    params1_true = extractRPSParams(A,1);
    params2_true = extractRPSParams(A,2);
    
    % Prepare the x tick labels
    x_txts = cell(1,length(simResults.pop_sizes));
    for n = 1:length(simResults.pop_sizes)
        x_txts{n} = num2str(simResults.pop_sizes(n));
    end
    x_txts{end+1} = '\infty';
    
    % Load in the data for each method
    for m = 1:N_methods
        
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
        % Load data
        load(this_filename, 'simResults');
               
        % Plot the root mean square error boxplot
        plotRMSEBoxplotOneType( struct('l2_mat',simResults.obsl2_mat,'l2_perfect',simResults.obsl2_perfect), x_txts, [problem_names{p},', ',method_txts{m}], ['FIGS_',problem_filenames{p},'_',method_filenames{m},'_observedRMSE'], save_figs, 1, true);
        plotRMSEBoxplotOneType( struct('l2_mat',simResults.obsl2_mat,'l2_perfect',simResults.obsl2_perfect), x_txts, [problem_names{p},', ',method_txts{m}], ['FIGS_',problem_filenames{p},'_',method_filenames{m},'_observedRMSE'], save_figs, 2, true);
        plotRMSEBoxplotOneType( struct('l2_mat',simResults.predl2_mat,'l2_perfect',simResults.predl2_perfect), x_txts, [problem_names{p},', ',method_txts{m}], ['FIGS_',problem_filenames{p},'_',method_filenames{m},'_forecastRMSE'], save_figs, 1, true);
        plotRMSEBoxplotOneType( struct('l2_mat',simResults.predl2_mat,'l2_perfect',simResults.predl2_perfect), x_txts, [problem_names{p},', ',method_txts{m}], ['FIGS_',problem_filenames{p},'_',method_filenames{m},'_forecastRMSE'], save_figs, 2, true);
                
        % Plot the parameter estimation boxplots
        plotParamsBoxplotOneType( setfield(simResults, 'params_true', params1_true), x_txts, '', ['FIG6b_',problem_filenames{p},'_',method_filenames{m},'_popsParams'], save_figs, 1, true );
        plotParamsBoxplotOneType( setfield(simResults, 'params_true', params2_true), x_txts, '', ['FIG6b_',problem_filenames{p},'_',method_filenames{m},'_popsParams'], save_figs, 2, true );
        
    end
    
end