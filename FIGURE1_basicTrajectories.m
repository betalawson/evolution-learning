function FIGURE1_basicTrajectories()

%%% SETUP

% Initialise random seed
rng(7);

% Create a function that parameterises the bi-allelic case in terms of only
% s and h, with base fitness of 1
allelePayoff = @(x, s, h) [1 + s, 1 + h * s; 1 + h*s, 1] * x;

% Specify the simulation parameters
selection_type = 1;
N_pop = 5000;
N_gen = 50;
x0 = [0.2;0.8];

% Specify the ranges over which to vary the two parameters
s_range = [0,1];
h_range = [-0.5, 1.5];

% Specify the number of curves to plot
N_plot = 11;

% Prepare a vector of generation numbers
r = 0:N_gen;

% Load the viridis colormap
load('extra_colormaps.mat','viridis');

% Determine how many evenly spaced values will fit
N_clr = size(viridis,1);
N_sep = floor( (N_clr-1) / (N_plot-1) );

% Section out these equally-spaced colours
traj_clrs = viridis(1:N_sep:N_clr,:);


%%% PLOTTING PREPARATION

% Number of slots occupied by each individual figure (increase to reduce
% gap between them)
N_slots = 11;
N_gap = 3;

% Prepare the figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
tiles_obj = tiledlayout(1,2*N_slots+2*N_gap);
tiles_obj.TileSpacing = 'loose';
tiles_obj.Padding = 'tight';

%%% SELECTION STRENGTH PLOT

% Move to next tile in tiledLayout
ax = nexttile([1 N_slots]);
ax.NextPlot = 'replacechildren';
hold on;

% Create the colorbar and make it wider
clrbar_obj = colorbar(ax);
title(clrbar_obj,'$s$','FontSize',28,'Interpreter','LaTeX');

% Set up parameters
s_vals = linspace(s_range(1), s_range(2), N_plot);
h = 0.8;

% Loop over varied parameter
for k = 1:N_plot

    % Generate a Wright-Fisher trajectory
    X = wrightFisher(N_pop, N_gen, x0, @(x) allelePayoff(x, s_vals(k), h), selection_type);

    % Plot this trajectory
    plot(r,X(1,:)', 'LineWidth', 3.5, 'Color', traj_clrs(k,:),'parent',ax);

end

% Carry out plot clean-up and add colorbar as "legend" for trajectories
xlabel('Generation Number, r');
ylabel('Proportion of Favoured Allele, x_1');
set(gca,'FontSize',24,'LineWidth',2);
title('h = 0.8,      s varying','FontSize',28);
axis square;
colormap(traj_clrs);
caxis(s_range);
xlim([0 N_gen]);
clrbar_obj.Position = [0.4, 0.3, 0.025, 0.45];

%%% EMPTY SPACE
nexttile([1 N_gap]);
set(gca,'Visible',false);


%%% DOMINANCE PLOT - COPY PASTE AFTER PRETTYING UP THE ABOVE

% Move to next tile in tiledLayout
ax = nexttile([1 N_slots]);
ax.NextPlot = 'replacechildren';
hold on;

% Create the colorbar and make it wider
clrbar_obj = colorbar(ax);
title(clrbar_obj,'$h$','FontSize',28,'Interpreter','LaTeX');

% Set up parameters
h_vals = linspace(h_range(1), h_range(2), N_plot); 
s = 0.2;

% Loop over varied parameter
for k = 1:N_plot

    % Generate a Wright-Fisher trajectory
    X = wrightFisher(N_pop, N_gen, x0, @(x) allelePayoff(x, s, h_vals(k)), selection_type);

    % Plot this trajectory
    plot(r,X(1,:)', 'LineWidth', 3.5, 'Color', traj_clrs(k,:),'parent',ax);

end

% Carry out plot clean-up and add colorbar as "legend" for trajectories
xlabel('Generation Number, r');
ylabel('Proportion of Favoured Allele, x_1');
set(gca,'FontSize',24,'LineWidth',2);
title('s = 0.2,      h varying','FontSize',28);
axis square;
colormap(traj_clrs);
caxis(h_range);
xlim([0 N_gen]);
clrbar_obj.Position = [0.9075, 0.3, 0.025, 0.45];