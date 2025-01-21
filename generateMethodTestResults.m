function generateMethodTestResults(regenerate)

% This function tests parameter estimation of standard replicator match,
% gradient fisher, and Wright-Fisher max likelihood with known N

%%% PARAMETER SETTINGS (ONLY USED IF REGENERATING THE DATA)

% Specify true values of the parameters
s_true_vals = linspace(0,0.25,6);             % Selective advantage
h_true = 0.8;                                 % Dominance in heterozygotes
X0 = [0.2;0.8];                               % Initial proportions

% Specify the data to which the parameters are being fit
N_pop = 1000;             % Population size for WF model
N_gen = 20;               % Number of generations in WF model (or t_end in replicator)
t_gen = 1;                % Time associated with each generation
selection_type = 2;       % Type I or type II selection

% Specify the number of generations to use in evaluating forecast error
N_forecast = 50;

% Specify the number of replicates for box plots
N_rep = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATA REGENERATION

% Assume data is not to be regenerated unless necessary, or requested
if nargin < 1
    regenerate = false;
end

% Equation learning options
ELoptions = addDefaultOptions(struct('symmetric_payoff',true,'library_orders',1),2);

% Set data filename
data_filename = 'DATA_methodTests.mat';

% Define the type II replicator RHS for the two allele case
fitness_fun = @(x,s,hs) [1 + s, 1 + hs; 1 + hs, 1] * x;

% Generate the data if not present or asked to regenerate
if regenerate || ~exist(data_filename, 'file')
    
    % Prepare data storage
    s_data = zeros(length(s_true_vals), 3, N_rep);
    hs_data = zeros(length(s_true_vals), 3, N_rep);
    obsl2_data = zeros(length(s_true_vals), 3, N_rep);
    predl2_data = zeros(length(s_true_vals), 3, N_rep);
    hs_true_vals = zeros(length(s_true_vals));
    
    % Loop over different selection strengths
    for ii = 1:length(s_true_vals)
        
        % Read out current selection
        s_true = s_true_vals(ii);
        
        % Prepare problem object
        problem = struct('fitness', @(X) fitness_fun(X,s_true,h_true*s_true), 'N_pop', N_pop, 'N_gen', N_gen, 't_gen', t_gen, 'X0', X0, 'selection_type', selection_type, 'N_feat', 2, 'obs_Nth', 1);
        
        % Set seed for consistency
        rng(7);
        
        % Prepare optimisation options
        optimiser_options = optimset('fminunc');
        optimiser_options.MaxFunEvals = 20000;
        optimiser_options.MaxIter = 10000;
        
        % Repeatedly generate Wright-Fisher data, and parameter estimate
        ls_s = zeros(N_rep,1);
        ls_hs = zeros(N_rep,1);
        gm_s = zeros(N_rep,1);
        gm_hs = zeros(N_rep,1);
        ll_s = zeros(N_rep,1);
        ll_hs = zeros(N_rep,1);
        
        parfor r = 1:N_rep
            
            % Generate the Wright-Fisher and replicator data
            traj = simulateTrajectories(problem);
            
            % Convert to non-dimensionalised time for equation learning
            WF_nonD = struct('t',traj.WF.t / problem.t_gen, 'X', traj.WF.X);
            
            % Learn fitness using least squares and estimate parameters
            ls_library = evolutionLearning(WF_nonD, selection_type, setfield(ELoptions,'regression_type','ls'));
            ls_params = extractBialleleParams( libraryToPayoff(ls_library, 2), selection_type);
            ls_s(r) = ls_params.s;
            ls_hs(r) = ls_params.hs;
            
            % Learn fitness using gradient matching
            gm_library = evolutionLearning(WF_nonD, selection_type, setfield(ELoptions,'regression_type','gm'));
            gm_params = extractBialleleParams( libraryToPayoff(gm_library, 2), selection_type);
            gm_s(r) = gm_params.s;
            gm_hs(r) = gm_params.hs;
            
            % Learn fitness using Wright-Fisher likelihood
            params = fminunc(@(params) -WFLikelihood2(WF_nonD, setfield(problem, 'fitness', @(X) fitness_fun(X, params(1), params(2)))) / 21 + 0.001 * sqrt( params(1)^2 + params(2)^2), [s_true;h_true*s_true], optimiser_options);
            ll_s(r) = params(1);
            ll_hs(r) = params(2);
            
            LL_fitness = @(x) fitness_fun(x,params(1),params(2));
            GM_fitness = @(X) evaluateLibrary(X, gm_library);
            LS_fitness = @(X) evaluateLibrary(X, ls_library);
            
            ls_obsl2(r) = calculateTrajectoryError(problem, LS_fitness, problem.X0);
            gm_obsl2(r) = calculateTrajectoryError(problem, GM_fitness, problem.X0);
            ll_obsl2(r) = calculateTrajectoryError(problem, LL_fitness, problem.X0);
            
            ls_predl2(r) = calculateTrajectoryError(setfield(problem,'N_gen',N_forecast), LS_fitness, problem.X0);
            gm_predl2(r) = calculateTrajectoryError(setfield(problem,'N_gen',N_forecast), GM_fitness, problem.X0);
            ll_predl2(r) = calculateTrajectoryError(setfield(problem,'N_gen',N_forecast), LL_fitness, problem.X0);
            
        end
        
        % Store all data ready for box plotting
        s_data(ii,1,:) = ls_s;
        s_data(ii,2,:) = gm_s;
        s_data(ii,3,:) = ll_s;
        hs_data(ii,1,:) = ls_hs;
        hs_data(ii,2,:) = gm_hs;
        hs_data(ii,3,:) = ll_hs;
        obsl2_data(ii,1,:) = ls_obsl2;
        obsl2_data(ii,2,:) = gm_obsl2;
        obsl2_data(ii,3,:) = ll_obsl2;
        predl2_data(ii,1,:) = ls_predl2;
        predl2_data(ii,2,:) = gm_predl2;
        predl2_data(ii,3,:) = ll_predl2;
        
        % Also store true value of hs
        hs_true_vals(ii) = h_true * s_true;
        
    end
    
    save(data_filename,'s_data','hs_data','obsl2_data','predl2_data','s_true_vals','hs_true_vals');
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l2 = calculateTrajectoryError(Pobj, F, X0list)
% This function generates replicator trajectories for the provided problem
% object, and for the provided object but with the specified fitness
% function F, and evaluates the average l2 error between them across all of
% the initial conditions provided in X0list

% Read out the number of ICs given (number of column vectors in X0list)
N_ICs = size(X0list,2);

% Loop over all trajectories and evaluate each's contribution to mean error
l2 = 0;
for k = 1:N_ICs
               
    % Create a problem object for this initial condition
    Phere = setfield(Pobj, 'X0', X0list(:,k));
    
    % Create trajectory for true selective dynamics
    true_traj = simulateTrajectories( Phere, 'rep' );
    % Create trajectory for learned dynamics
    learned_traj = simulateTrajectories( setfield(Phere, 'fitness', F), 'rep' );
                
    % Calculate error - if integration failed despite adaptive tolerance
    % specification, store a NaN
    if isequal(size(true_traj.rep.X(:)),size(learned_traj.rep.X(:))) 
        l2 = l2 + rms( true_traj.rep.X(:) - learned_traj.rep.X(:) ) / N_ICs;
    else
        l2 = NaN;
    end
                
end

end
