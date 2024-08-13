function FIGURE_basicTrajectories()

%%% SETUP

% Create a function that parameterises the bi-allelic case in terms of only
% s and h, with base fitness of 1
allelePayoff = @(x, s, h) [1 + s, 1 + h * s; 1 + h*s, 1] * x;

% Specify the simulation parameters
selection_type = 1;
N_pop = 10000;
N_gen = 50;
x0 = [0.2;0.8];

% Specify the ranges over which to vary the two parameters
s_range = [0,1];
h_range = [-0.5, 1.5];

% Specify the number of curves to plot
N_plot = 11;

% Load the viridis colormap
load('extra_colormaps.mat','viridis');

% Determine how many evenly spaced values will fit
N_clr = size(viridis,1);
N_sep = floor( (N_clr-1) / N_plot);

% Section out these equally-spaced colours
traj_clrs = viridis(1:N_sep:N_clr,:);


%%% SELECTION STRENGTH PLOT - PRETTY UP LATER

% Set up parameters
s_vals = linspace(s_range(1), s_range(2), N_plot); 
h = 0.8;

% Prepare figure
figure; hold on;

% Loop over varied parameter
for k = 1:N_plot

    % Generate a Wright-Fisher trajectory
    X = wrightFisher(N_pop, N_gen, x0, @(x) allelePayoff(x, s_vals(k), h), selection_type);

    % Plot this trajectory
    plot(X(1,:)', 'LineWidth', 2, 'Color', traj_clrs(k,:));

end

% Carry out plot clean-up and add colorbar as "legend" for trajectories
xlabel('Number of Generations');
ylabel('Proportion of Favoured Allele');
set(gca,'FontSize',20);
colormap(traj_clrs);
caxis(s_range);
clrbar_obj = colorbar;
title(clrbar_obj,'$s$','FontSize',24,'Interpreter','LaTeX');


%%% DOMINANCE PLOT - COPY PASTE AFTER PRETTYING UP THE ABOVE

% Set up parameters
h_vals = linspace(h_range(1), h_range(2), N_plot); 
s = 0.2;

% Prepare figure
figure; hold on;

% Loop over varied parameter
for k = 1:N_plot

    % Generate a Wright-Fisher trajectory
    X = wrightFisher(N_pop, N_gen, x0, @(x) allelePayoff(x, s, h_vals(k)), selection_type);

    % Plot this trajectory
    plot(X(1,:)', 'LineWidth', 2, 'Color', traj_clrs(k,:));

end

% Carry out plot clean-up and add colorbar as "legend" for trajectories
xlabel('Number of Generations');
ylabel('Proportion of Favoured Allele');
set(gca,'FontSize',20);
colormap(traj_clrs);
caxis(h_range);
clrbar_obj = colorbar;
title(clrbar_obj,'$h$','FontSize',24,'Interpreter','LaTeX');