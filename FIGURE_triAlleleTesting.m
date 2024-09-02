function FIGURE_triAlleleTesting(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying selection strength and
% dominance for a bi-allelic inheritance scenario

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

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
load('problems_inheritance.mat', 'allele_standard1', 'allele_transient1', 'allele_persistent1', 'allele_standard2', 'allele_transient2', 'allele_persistent2');

% Create a cell array of cell arrays that contains the type I and type II
% versions of each problem to be looped over
problems_list = { {allele_standard1, allele_standard2}, {allele_transient1, allele_transient2}, {allele_persistent1, allele_persistent2} };

% Specify the filename "suffixes" for each problem in the master list
problem_filenames = {'std', 'tran', 'pers'};
problem_names = {'Standard', 'Transient', 'Persistent'};

% Loop over and run each problem from the problem list separately
for P = 1:length(problems_list)
    
    % Load in current problems (type I and type II) from master list
    problems = problems_list{P};
        
    %%% PASTE OTHER METHODS HERE IF THEY ARE ALSO TO BE RUN
    
    % Using method: Symmetric payoff matrix
    filename = ['DATA_IN',problem_filenames{P},'_Npop_symm'];
    % Regenerate if requested or data not present
    if ~exist([filename,'.mat'],'file') || regenerate
        ELoptions.symmetric_payoff = true;
        ELoptions.library_orders = [1];
        simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
    else
        load([filename,'.mat'],'simResults');
    end    
    plotOverPopSizes(simResults, problem_names{P});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename)
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
s1_mat = zeros(N_pops, 2, N_rep);
s2_mat = zeros(N_pops, 2, N_rep);
h12_mat = zeros(N_pops, 2, N_rep);
h13_mat = zeros(N_pops, 2, N_rep);
h23_mat = zeros(N_pops, 2, N_rep);
l2_mat = zeros(N_pops, 2, N_rep);

for k = 1:N_pops
   
    % Loop over the two types of selection and generate results for each
    for type = [1,2]
    
        % Overwrite the population size for this problem
        problem = problems{type};
        problem.N_pop = pop_sizes(k);
    
        % Repeatedly generate new Wright-Fisher data and learn s and h
        for r = 1:N_rep
        
            % Generate the Wright-Fisher data
            traj = generateTrajectories(problem);
        
            % Convert to non-dimensionalised time for equation learning
            WF_nonD = struct('t',traj.WF.t / problem.t_gen, 'X', traj.WF.X);
        
            % Build additional trajectories for alternate initial
            % conditions as extra data
%              IClist = [ [0.1;0.4;0.5], [0.2; 0.7; 0.1], [0.4; 0.3; 0.3], [0.4; 0.5; 0.1], [0.4; 0.1; 0.5], [0.6; 0.3; 0.1], [0.6; 0.1; 0.3], [0.1; 0.8; 0.1], [0.1; 0.6; 0.3], [0.1; 0.3; 0.6] ];
%              traj_data = cell(1,size(IClist,2)+1);
%              traj_data{1} = WF_nonD;
%              for kk = 1:size(IClist,2)
%                  ICproblem = problem;
%                  ICproblem.X0 = IClist(:,kk);
%                  ICtraj = generateTrajectories(ICproblem);
%                  ICWF_nonD = struct('t',ICtraj.WF.t / problem.t_gen, 'X', ICtraj.WF.X);
%                  traj_data{kk+1} = ICWF_nonD;
%              end
            
            % Perform equation learning
            [fitfun, Flib, K, texts] = evolutionLearning(WF_nonD, type, ELoptions);
        
            % Extract learned payoff matrix and normalise using a_33 = 1
            payoff = fitnessToPayoff(K,texts,3);
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
                b = [ payoff(1,2) - payoff(2,1) ; payoff(1,3) - payoff(3,1) ; payoff(2,3) - payoff(3,2) ];
                % Find minimum-norm x that gets as close as possible to
                % satisfying symmetry checks
                z = lsqminnorm(C, b);

                % Adjust the payoff matrix using these shifts
                payoff = payoff + z';

                % Now adjust the payoff matrix so the bottom-right element
                % is unity to match the base matrix definition
                payoff = payoff - payoff(3,3) + 1;

            else
                
                % Type II selection only allows constant scaling, so we
                % cannot do anything to push towards symmetry, only
                % normalise using bottom-right element
                payoff = payoff / payoff(3,3);

            end

            % Get equivalent symmeterised off-diagonal elements
            a12 = 0.5 * ( payoff(1,2) + payoff(2,1) );
            a13 = 0.5 * ( payoff(1,3) + payoff(3,1) );
            a23 = 0.5 * ( payoff(2,3) + payoff(3,2) );

            % Extract the values of s and h and store them 
            s1 = payoff(1,1) - 1;
            s2 = payoff(2,2) - 1;
            h12 = ( ( a12 - 1 ) - s2 ) / (s1 - s2);
            h13 = ( a13 - 1 ) / s1;
            h23 = ( a23 - 1 ) / s2;
            
            % Store these in their respective matrices
            s1_mat(k,type,r) = s1;
            s2_mat(k,type,r) = s2;
            h12_mat(k,type,r) = h12;
            h13_mat(k,type,r) = h13;
            h23_mat(k,type,r) = h23;
            
            % Simulate the learned dynamics and evaluate the L2 error
            % between the learned dynamics and the actual replicator
            learned_problem = problem;
            learned_problem.fitness = @(X) payoff*X;
            learned_traj = generateTrajectories(learned_problem, 'rep');
            l2_mat(k,type,r) = rms(traj.rep.X(:) - learned_traj.rep.X(:));      

        end 
    end    
    
end

% Now generate "perfect" replicator data with known derivative and estimate
% using it for reference
s1_perfect = zeros(1,2);
s2_perfect = zeros(1,2);
h12_perfect = zeros(1,2);
h13_perfect = zeros(1,2);
h23_perfect = zeros(1,2);
l2_perfect = zeros(1,2);

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
    
    % Extract learned payoff matrix and normalise using a_33 = 1
    payoff = fitnessToPayoff(K,texts,3);
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
        b = [ payoff(1,2) - payoff(2,1) ; payoff(1,3) - payoff(3,1) ; payoff(2,3) - payoff(3,2) ];
        % Find minimum-norm x that gets as close as possible to
        % satisfying symmetry checks
        z = lsqminnorm(C, b);

        % Adjust the payoff matrix using these shifts
        payoff = payoff + z';

        % Now adjust the payoff matrix so the bottom-right element
        % is unity to match the base matrix definition
        payoff = payoff - payoff(3,3) + 1;

    else

        % Type II selection only allows constant scaling, so we
        % cannot do anything to push towards symmetry, only
        % normalise using bottom-right element
        payoff = payoff / payoff(3,3);

    end

    % Get equivalent symmeterised off-diagonal elements
    a12 = 0.5 * ( payoff(1,2) + payoff(2,1) );
    a13 = 0.5 * ( payoff(1,3) + payoff(3,1) );
    a23 = 0.5 * ( payoff(2,3) + payoff(3,2) );

    % Extract the values of s and h and store them
    s1 = payoff(1,1) - 1;
    s2 = payoff(2,2) - 1;
    h12 = ( ( a12 - 1 ) - s2 ) / (s1 - s2);
    h13 = ( a13 - 1 ) / s1;
    h23 = ( a23 - 1 ) / s2;
    s1_perfect(type) = s1;
    s2_perfect(type) = s2;
    h12_perfect(type) = h12;
    h13_perfect(type) = h13;
    h23_perfect(type) = h23;
    
    % Find the L2 error to the true selective trajectory
    learned_problem = problem;
    learned_problem.fitness = @(X) payoff*X;
    learned_traj = generateTrajectories(learned_problem);
    l2_perfect(type) = rms(traj.rep.X(:) - learned_traj.rep.X(:));

end

% Extract the true values of s and h from the problem
payoff = problems{1}.fitness(eye(length(problem.X0)));
% Matrix should already be in correct form, but do all the equivalent
% normalisations anyway as per the above
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
    b = [ payoff(1,2) - payoff(2,1) ; payoff(1,3) - payoff(3,1) ; payoff(2,3) - payoff(3,2) ];
    % Find minimum-norm x that gets as close as possible to
    % satisfying symmetry checks
    z = lsqminnorm(C, b);

    % Adjust the payoff matrix using these shifts
    payoff = payoff + z';

    % Now adjust the payoff matrix so the bottom-right element
    % is unity to match the base matrix definition
    payoff = payoff - payoff(3,3) + 1;

else

    % Type II selection only allows constant scaling, so we
    % cannot do anything to push towards symmetry, only
    % normalise using bottom-right element
    payoff = payoff / payoff(3,3);

end
    
% Get equivalent symmeterised off-diagonal elements
a12 = 0.5 * ( payoff(1,2) + payoff(2,1) );
a13 = 0.5 * ( payoff(1,3) + payoff(3,1) );
a23 = 0.5 * ( payoff(2,3) + payoff(3,2) );

% Extract the values of s and h and store them
s1 = payoff(1,1) - 1;
s2 = payoff(2,2) - 1;
h12 = ( ( a12 - 1 ) - s2 ) / (s1 - s2);
h13 = ( a13 - 1 ) / s1;
h23 = ( a23 - 1 ) / s2;
s1_true = s1;
s2_true = s2;
h12_true = h12;
h13_true = h13;
h23_true = h23;

% Store all data in a single results object and save it
simResults = struct('s1_mat',s1_mat,'s2_mat',s2_mat,'h12_mat',h12_mat,'h13_mat',h13_mat,'h23_mat',h23_mat,'l2_mat',l2_mat,'s1_perfect',s1_perfect,'s2_perfect',s2_perfect,'h12_perfect',h12_perfect,'h13_perfect',h13_perfect,'h23_perfect',h23_perfect,'l2_perfect',l2_perfect,'s1_true',s1_true,'s2_true',s2_true,'h12_true',h12_true,'h13_true',h13_true,'h23_true',h23_true,'pop_sizes',pop_sizes);
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
s1_mat = simResults.s1_mat;
s2_mat = simResults.s2_mat;
h12_mat = simResults.h12_mat;
h13_mat = simResults.h13_mat;
h23_mat = simResults.h23_mat;
l2_mat = log(simResults.l2_mat);
s1_perfect = simResults.s1_perfect;
s2_perfect = simResults.s2_perfect;
h12_perfect = simResults.h12_perfect;
h13_perfect = simResults.h13_perfect;
h23_perfect = simResults.h23_perfect;
l2_perfect = log(simResults.l2_perfect);
s1_true = simResults.s1_true;
s2_true = simResults.s2_true;
h12_true = simResults.h12_true;
h13_true = simResults.h13_true;
h23_true = simResults.h23_true;
pop_sizes = simResults.pop_sizes;

% Gather data together into cell arrays that enable the below for loop
p_mats = {s1_mat, s2_mat, h12_mat, h13_mat, h23_mat};
p_trues = {s1_true, s2_true, h12_true, h13_true, h23_true};
p_perfects = {s1_perfect, s2_perfect, h12_perfect, h13_perfect, h23_perfect};

% Specify the parameter names (letters and descriptions)
p_names = {'s_1', 's_2', 'h_{12}', 'h_{13}', 'h_{23}'};
p_descs = {'Selective Advantage of Allele 1', 'Selective Advantage of Allele 2', 'Dominance of Allele 1 over 2', 'Dominance of Allele 1 over 3', 'Dominance of Allele 2 over 3'}; 

% Prepare plot
N_pops = length(pop_sizes);

% Prepare the labels for each population size
x_txts = cell(1,N_pops);
for k = 1:N_pops
    x_txts{k} = ['N = ',num2str(pop_sizes(k))];
end
x_txts{N_pops+1} = 'Perfect Data';


% Loop over each bit of data in the cell arrays
for pnum = 1:length(p_mats)
    
    % Extract this parameter's data
    p_mat = p_mats{pnum};
    p_true = p_trues{pnum};
    p_perfect = p_perfects{pnum};

    % Prepare this figure
    figure; hold on;
    % Plot the true value first
    plot([-1 N_pops+2],[p_true p_true], 'k', 'LineWidth', 1.5);
    % Plot the boxplot showing main data
    boxplot_obj = boxplot2(p_mat);
    % Add colours to the plot
    for ii = 1:2
        structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
    end
    % Add to the plot the parameter learned for perfect data
    hold on;
    plot(N_pops+0.9,p_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
    plot(N_pops+1.1,p_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
    % Clean up plot
    axis([0 N_pops+1.2 ymin ymax]);
    xticks(1:N_pops+1);
    xticklabels(x_txts);
    set(gca,'FontSize',20);
    ylabel(p_names{pnum},'Fontsize',24);
    title(['EstiMATion of ',p_descs{pnum},' --- ',title_txt],'FontSize',20);

end


%%% SEPARATE FIGURE FOR THE ERROR BETWEEN LEARNED DYNAMICS AND REPLICATOR

% Prepare the figure
figure; hold on;
% Plot the boxplot showing main data
boxplot_obj = boxplot2(l2_mat);
% Add colours to the plot
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', type_clrs(ii,:), 'markeredgecolor', type_clrs(ii,:), 'LineWidth', 2), boxplot_obj);
end
hold on;
plot(N_pops+0.9,l2_perfect(1),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(1,:));
plot(N_pops+1.1,l2_perfect(2),'.','MarkerSize',40, 'MarkerEdgeColor',type_clrs(2,:));
% Clean up plot
axis([0 N_pops+1.2  -0.1+min([l2_mat(:);l2_perfect(:)]) 0.1+max([l2_mat(:); l2_perfect(:)]) ]);
xticks(1:N_pops+1);
xticklabels(x_txts);
set(gca,'FontSize',20);
ylabel('Root Mean Square Error','Fontsize',24);
title(['L_2 error to True Selection Trajectory --- ',title_txt],'FontSize',20);


% THE BELOW MAY BE COPY-PASTED INTO THE SECTION SPECIFIED ABOVE (WITH MINOR
% MODIFICATIONS OF TITLE TEXT ON FIGURES ETC) TO TEST OTHER METHODS THAT
% DON'T GUARANTEE A SYMMETRIC PAYOFF MATRIX

% %%% RUNNING PROBLEMS UNDER DIFFERENT SETTINGS
% 
% % General linear equation learning
% title_txt = 'Up to First Order';
% filename = ['DATA_IN',problem_filenames{P},'_Npop_orders01'];
% % Regenerate if requested or data not present
% if ~exist([filename,'.mat'],'file') || regenerate
%     ELoptions.symmetric_payoff = false;
%     ELoptions.library_orders = [0,1];
%     simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
% else
%     load([filename,'.mat'],'simResults');
% end    
% plotOverPopSizes(simResults, title_txt);
% 
% % Strictly payoff matrix elements
% title_txt = 'Strictly First Order';
% filename = ['DATA_IN',problem_filenames{P},'_Npop_orders1'];
% % Regenerate if requested or data not present
% if ~exist([filename,'.mat'],'file') || regenerate
%     ELoptions.symmetric_payoff = false;
%     ELoptions.library_orders = [1];
%     simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
% else
%     load([filename,'.mat'],'simResults');
% end    
% plotOverPopSizes(simResults, title_txt);
% 
% % Quadratic fitness function
% title_txt = 'Up to Second Order';
% filename = ['DATA_IN',problem_filenames{P},'_Npop_orders012'];
% % Regenerate if requested or data not present
% if ~exist([filename,'.mat'],'file') || regenerate
%     ELoptions.symmetric_payoff = false;
%     ELoptions.library_orders = [0,1,2];
%     simResults = runNpopSims(problems, N_rep, pop_sizes, ELoptions, filename);
% else
%     load([filename,'.mat'],'simResults');
% end    
% plotOverPopSizes(simResults, title_txt);
