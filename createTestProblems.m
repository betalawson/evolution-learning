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
N_gen = 200;
% Time associated with each generation
t_gen = 10;
% Observation frequency for WF model data
obs_Nth = 10;

% Allele selection strength and dominance
s = 0.2;
h = 0.8;
fitness = fitnessAlleles([s;0],[0,0;h,0]);

% Create problem objects
allele_basic1 = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth,'fitness',fitness,'selection_type',1);
allele_basic2 = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth,'fitness',fitness,'selection_type',2);

% Save problem objects
save('problems_basic.mat','allele_basic1','allele_basic2');


%%%%% ALLELE INHERITANCE PROBLEMS %%%%%

%%% CONSISTENT SETTINGS FOR ALL PROBLEMS

% Number of features
N_feat = 3;
% Initial condition
X0 = [0.1; 0.1; 0.8];
% Population size
N_pop = 2000;
% Number of generations
N_gen = 200;
% Time associated with each generation
t_gen = 10;
% Observation frequency for WF model data
obs_Nth = 10;

% Prepare the base problem used to construct the specific scenarios
allele_problem = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth);


%%% STANDARD TRI-ALLELIC SELECTION

% Feature selection strengths
s_1 = 1;
s_2 = 0.5;
% Heterozygote dominance values
h_12 = 0.5;
h_13 = 1;
h_23 = 1;
% Generate fitness function
fitness = fitnessAlleles([s_1,s_2,0], [0,0,0; h_12, 0, 0; h_13, h_23, 0]);

% Create type I and type II selection versions of this problem
allele_standard = allele_problem;
allele_standard.fitness = fitness;
allele_standard1 = allele_standard;
allele_standard1.selection_type = 1;
allele_standard2 = allele_standard;
allele_standard2.selection_type = 2;


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

% Create type I and type II selection versions of this problem
allele_transient = allele_problem;
allele_transient.fitness = fitness;
allele_transient1 = allele_transient;
allele_transient1.selection_type = 1;
allele_transient2 = allele_transient;
allele_transient2.selection_type = 2;


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

% Create type I and type II selection versions of this problem
allele_persistent = allele_problem;
allele_persistent.fitness = fitness;
allele_persistent1 = allele_persistent;
allele_persistent1.selection_type = 1;
allele_persistent2 = allele_persistent;
allele_persistent2.selection_type = 2;


%%% SAVE THE PROBLEMS FOR PLOTTING OR FOR APPLYING EQUATION LEARNING
save('problems_inheritance.mat', 'allele_standard1', 'allele_transient1', 'allele_persistent1', 'allele_standard2', 'allele_transient2', 'allele_persistent2');

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
N_gen = 200;
% Time associated with each generation
t_gen = 10;
% Observation frequency for WF model data
obs_Nth = 1;

% Prepare the base problem used to construct the specific scenarios
RPS_problem = struct('N_feat',N_feat,'X0',X0,'N_pop',N_pop,'N_gen',N_gen,'t_gen',t_gen,'obs_Nth',obs_Nth);


%%% BALANCED ROCK PAPER SCISSORS

% Initialise problem
RPS_balanced = RPS_problem;
% Prepare balanced version of RPS dynamics - b1b2b3 = c1c2c3
f_base = 2;
win_strength = [0.3, 0.8, 0.2];
loss_strength = [0.3, 0.8, 0.2];
RPS_balanced.fitness = fitnessRPS(N_feat, f_base, win_strength, loss_strength);
% Create type I and type II selection versions of this problem
RPS_balanced1 = RPS_balanced;
RPS_balanced1.selection_type = 1;
RPS_balanced2 = RPS_balanced;
RPS_balanced2.selection_type = 2;


%%% ATTRACTING ROCK PAPER SCISSORS

% Initialise problem
RPS_attract = RPS_problem;
% Prepare attracting version of RPS dynamics - b1b2b3 > c1c2c3
f_base = 2;
win_strength = [0.2, 1.8, 0.2];
loss_strength = [0.3, 0.8, 0.2];
RPS_attract.fitness = fitnessRPS(N_feat, f_base, win_strength, loss_strength);
% Create type I and type II selection versions of this problem
RPS_attract1 = RPS_attract;
RPS_attract1.selection_type = 1;
RPS_attract2 = RPS_attract;
RPS_attract2.selection_type = 2;


%%% REPELLING ROCK PAPER SCISSORS

% Initialise problem
RPS_repel = RPS_problem;
% Prepare attracting version of RPS dynamics - b1b2b3 > c1c2c3
f_base = 2;
win_strength = [0.3, 0.3, 0.2];
loss_strength = [0.3, 0.5, 0.2];
RPS_repel.fitness = fitnessRPS(N_feat, f_base, win_strength, loss_strength);
% Create type I and type II selection versions of this problem
RPS_repel1 = RPS_repel;
RPS_repel1.selection_type = 1;
RPS_repel2 = RPS_repel;
RPS_repel2.selection_type = 2;


%%% SAVE THE PROBLEMS FOR PLOTTING OR FOR APPLYING EQUATION LEARNING
save('problems_RPS.mat', 'RPS_balanced1', 'RPS_attract1', 'RPS_repel1', 'RPS_balanced2', 'RPS_attract2', 'RPS_repel2');

end