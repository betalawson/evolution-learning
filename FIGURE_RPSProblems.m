function FIGURE_RPSProblems()
% This function loads in the trajectories for the rock-paper-scissors
% problems and visualises them, for display in the paper. Note that the
% problems themselves will have been already created using the function
% createTestProblems.m and those problem objects will include properties
% such as the type of selection and the population size.

% Load in the problems for allele inheritance
load('problems_RPS.mat', 'RPS_balanced1', 'RPS_attract1', 'RPS_repel1', 'RPS_balanced2', 'RPS_attract2', 'RPS_repel2');
problems = {RPS_balanced2, RPS_attract2, RPS_repel2};

% Specify the titles for the figures
titles = {'Balanced ($b_1b_2b_3 = c_1c_2c_3$)','Attracting ($b_1b_2b_3 > c_1c_2c_3$)','Repelling ($b_1b_2b_3 < c_1c_2c_3$)'};

% Set random seed for consistent Wright-Fisher trajectories
rng(7);
% Ensure we observe every generation in WF data for visualisation
for k = 1:length(problems)
    problems{k}.obs_Nth = 1;
end
% Generate the trajectories of each type I problem (WF and replicator)
RPS_trajectories = generateTrajectories( problems );

% Visualise these scenarios
for k = 1:length(RPS_trajectories)
    
    % Grab current trajectory from loop
    traj = RPS_trajectories{k};
    
    % First plot gets a legend, others do not
    if k == 1
        plotTriangleTrajectory(traj.WF.X, traj.rep.X, titles{k}, true);
    else
        plotTriangleTrajectory(traj.WF.X, traj.rep.X, titles{k}, false);
    end

end