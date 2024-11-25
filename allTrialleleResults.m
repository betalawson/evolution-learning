function allTrialleleResults(regenerate)
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

% Specify how many generations to simulate for prediction error checking
sim_gen = 50;

% Dataset name
data_prefix = 'DATA_TRIALL_';

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

% Load in the tri-allelic selection problems
load('problems_inheritance.mat', 'allele_standard', 'allele_transient', 'allele_persistence');

% Create a cell array of cell arrays that contains the problem list
problems_list = { allele_standard, allele_transient, allele_persistence };

% Specify the filename "suffixes" for each problem in the master list
problem_filenames = {'std', 'tran', 'pers'};
problem_names = {'Standard', 'Transient', 'Persistent'};

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2], [1] };
symmetric_flags = { false, false, false, true };

% Set up the filenames where data will be saved/loaded from
method_filenames = {'orders01','orders1','orders012','symm'};

% Count provided number of methods and number of problems
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
        ELoptions.symmetric_payoff = symmetric_flags{m};
    
        % Generate data for this problem if data not present or requested to
        if ~exist(this_filename,'file') || regenerate
            simResults = runSims( problems_list{p}, ELoptions, @(library,type) extractTrialleleParams(libraryToPayoff(library,3), type), pop_sizes, N_rep, sim_gen );
            save(this_filename, 'simResults', 'ELoptions');
        end
   
    end
end

% Safety pause
pause(1);


%%% PLOT FIGURE 4A: TRAJECTORY FIT QUALITY

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
    x_txts = {'p \leq 1', 'p = 1', 'p \leq 2', 'Symm'};
    title_txt = [problem_names{p},' N = ',num2str(N_pop)];
    % Plot the error boxplots
    plotRMSEBoxplot( struct('l2_mat', obsl2_mat), x_txts, 'Library', title_txt, ['FIG4a_',problem_filenames{p},'_observedRMSE'], save_figs, false, false);
    plotRMSEBoxplot( struct('l2_mat', predl2_mat), x_txts, 'Library', title_txt, ['FIG4a_',problem_filenames{p},'_forecastRMSE'], save_figs, false, false);
    
end


%%% PLOT FIGURE 4B: RECOVERY OF PARAMETERS FOR DIFFERENT POPULATION SIZES

% Loop over all problems
for p = 1:N_problems
    
    % Load in the data for the symmetric case
    load([data_prefix, problem_filenames{p},'_symm.mat'], 'simResults');

    % Prepare the x tick labels
    x_txts = cell(1,length(simResults.pop_sizes));
    for n = 1:length(simResults.pop_sizes)
        x_txts{n} = num2str(simResults.pop_sizes(n));
    end
    x_txts{end+1} = '\infty';

    % Get the true values of the parameters and add these to the results struct
    A = problems_list{p}.fitness(eye(3));
    simResults.params_true = extractTrialleleParams(A,1);

    plotParamsBoxplot( simResults, x_txts, '$N$', '', ['FIG4b_',problem_filenames{p},'popsParams'], save_figs, false, true );

end


%%% EXTRA PLOTS OF ERROR FOR DIFFERENT METHODS

% Provide text for methods
method_txts = {'$p \leq 1$', '$p = 1$', '$p \leq 2$', 'Symm.'};

% Loop over all problems
for p = 1:N_problems
    
    % Load in the data for each method
    for m = 1:N_methods
        
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
        % Load data
        load(this_filename, 'simResults');
        
        % Plot the root mean square error boxplot
        plotRMSEBoxplot( struct('l2_mat',simResults.obsl2_mat,'l2_perfect',simResults.obsl2_perfect), x_txts, 'Library', [problem_names{p},', ',method_txts{m}], ['FIGS_',problem_filenames{p},'_',method_filenames{m},'_observedRMSE'], save_figs, false, true);
        plotRMSEBoxplot( struct('l2_mat',simResults.predl2_mat,'l2_perfect',simResults.predl2_perfect), x_txts, 'Library', [problem_names{p},', ',method_txts{m}], ['FIGS_',problem_filenames{p},'_',method_filenames{m},'_forecastRMSE'], save_figs, false, true);
        
    end
    
end