function FIGURE_alleleProblems()
% This function loads in the trajectories for the allele inheritance
% problems and visualises them, for display in the paper. Note that the
% problems themselves will have been already created using the function
% createTestProblems.m and those problem objects will include properties
% such as the type of selection and the population size.

% Load in the problems for allele inheritance
load('problems_inheritance.mat', 'allele_standard1', 'allele_transient1', 'allele_persistent1');
problems = {allele_standard1, allele_transient1, allele_persistent1};

% Specify the titles for the figures
titles = {'Inheritance --- Standard','Inheritance --- Transient','Inheritance --- Persistence'};

% Set random seed for consistent Wright-Fisher trajectories
rng(7);
% Ensure we observe every generation in WF data for visualisation
for k = 1:length(problems)
    problems{k}.obs_Nth = 1;
end
% Generate the trajectories of each type I problem (WF and replicator)
allele_trajectories = generateTrajectories( problems );

% Visualise these scenarios
for k = 1:length(allele_trajectories)
    
    % Grab current trajectory from loop
    traj = allele_trajectories{k};
    
    % First plot gets a legend, others do not
    if k == 1
        plotTriangleTrajectory(traj.WF.X, traj.rep.X, titles{k}, true);
    else
        plotTriangleTrajectory(traj.WF.X, traj.rep.X, titles{k}, false);
    end

end