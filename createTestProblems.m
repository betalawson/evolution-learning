function createTestProblems()
% This function creates the test problems that are used in the paper. This
% function also visualises the problems, so they can be displayed in the
% paper.


%%%%% BASIC IDENTIFICATION OF SELECTION STRENGTH AND DOMINANCE

% Number of features
N_feat = 2;
% Initial condition
X0 = [0.2; 0.8];
% Population size
N_pop = 2000;
% Number of generations
N_gen = 20;
% Time associated with each generation
t_gen = 10;
% Observation frequency for WF model data
obs_Nth = 1;
% Selection type
selection_type = 1;

% Allele selection strength and dominance
s = 0.2;
h = 0.8;
fitness = fitnessAlleles([s;0],[0,0;h,0], 1);

% Create problem objects
allele_basic = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth,'fitness',fitness,'selection_type',selection_type);

% Save problem objects
save('problems_basic.mat','allele_basic');


%%%%% ALLELE INHERITANCE PROBLEMS %%%%%

%%% CONSISTENT SETTINGS FOR ALL PROBLEMS

% Number of features
N_feat = 3;
% Initial condition
X0 = [0.1; 0.1; 0.8];
% Population size
N_pop = 2000;
% Number of generations
N_gen = 20;
% Time associated with each generation
t_gen = 10;
% Observation frequency for WF model data
obs_Nth = 1;
% Selection type
selection_type = 1;

% Prepare the base problem used to construct the specific scenarios
allele_problem = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth,'selection_type',selection_type);


%%% STANDARD TRI-ALLELIC SELECTION

% Feature selection strengths
s_1 = 1;
s_2 = 0.5;
% Heterozygote dominance values
h_12 = 0.5;
h_13 = 1;
h_23 = 1;
% Generate fitness function
fitness = fitnessAlleles([s_1,s_2,0], [0, 0, 0; h_12, 0, 0; h_13, h_23, 0]);

% Add fitness to base problem to create this problem
allele_standard = setfield(allele_problem,'fitness',fitness);



%%% SELECTION EXHIBITING TRANSIENT SURGE 

% Feature selection strengths
s_1 = 1;
s_2 = 0.9;
% Heterozygote dominance values
h_12 = 0.5;
h_13 = 0.5;
h_23 = 1;
% Generate fitness function
fitness = fitnessAlleles([s_1,s_2,0], [0,0,0; h_12, 0, 0; h_13, h_23, 0]);

% Add fitness to base problem to create this problem
allele_transient = setfield(allele_problem,'fitness',fitness);


%%% SELECTION THAT PRODUCES PERSISTENCE OF ALL ALLELES

% Feature selection strengths
s_1 = 0.2;
s_2 = 0.15;
% Heterozygote dominance values
h_12 = 1.5;
h_13 = 1.25;
h_23 = 1.25;
% Generate fitness function
fitness = fitnessAlleles([s_1,s_2,0], [0,0,0; h_12, 0, 0; h_13, h_23, 0]);

% Add fitness to base problem to create this problem
allele_persistence = setfield(allele_problem,'fitness',fitness);


%%% SAVE THE PROBLEMS FOR PLOTTING OR FOR APPLYING EQUATION LEARNING
save('problems_inheritance.mat', 'allele_standard', 'allele_transient', 'allele_persistence');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% ROCK PAPER SCISSORS PROBLEMS %%%%%

%%% CONSISTENT SETTINGS FOR ALL PROBLEMS

% Number of features
N_feat = 3;
% Initial condition
X0 = [0.2; 0.5; 0.3];
% Population size
N_pop = 2000;
% Number of generations
N_gen = 100;
% Time associated with each generation
t_gen = 10;
% Observation frequency for WF model data
obs_Nth = 1;
% Selection type
selection_type = 1;

% Prepare the base problem used to construct the specific scenarios
RPS_problem = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth,'selection_type',selection_type);


%%% BALANCED ROCK PAPER SCISSORS

% Prepare balanced version of RPS dynamics: b1b2b3 = c1c2c3
f_base = 2;
win_strength = [0.3, 0.8, 0.2];
loss_strength = [0.3, 0.8, 0.2];
fitness = fitnessRPS(N_feat, f_base, win_strength, loss_strength);
% Create the problem with this fitness function
RPS_balanced = setfield(RPS_problem,'fitness',fitness);


%%% ATTRACTING ROCK PAPER SCISSORS

% Prepare attracting version of RPS dynamics - b1b2b3 > c1c2c3
f_base = 2;
win_strength = [0.3, 0.8, 0.2];
loss_strength = [0.3, 0.1, 0.2];
fitness = fitnessRPS(N_feat, f_base, win_strength, loss_strength);
% Create the problem with this fitness function
RPS_attract = setfield(RPS_problem,'fitness',fitness);


%%% REPELLING ROCK PAPER SCISSORS

% Prepare attracting version of RPS dynamics - b1b2b3 < c1c2c3
f_base = 2;
win_strength = [0.3, 0.3, 0.2];
loss_strength = [0.3, 0.5, 0.2];
fitness = fitnessRPS(N_feat, f_base, win_strength, loss_strength);
% Create the problem with this fitness function
RPS_repel = setfield(RPS_problem,'fitness',fitness);


%%% SAVE THE PROBLEMS FOR PLOTTING OR FOR APPLYING EQUATION LEARNING
save('problems_RPS.mat', 'RPS_balanced', 'RPS_attract', 'RPS_repel');

end