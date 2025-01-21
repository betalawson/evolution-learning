function FIGURES2_biallele_lossLandscapes(regenerate)
% This function plots the loss surface to be optimised in fitting the
% parameters of a basic bi-allelic inheritance problem to synthetic data.
% This is compared with the corresponding loss surface for gradient
% matching the corresponding derivative data.


%%% PARAMETER SETTINGS (ONLY USED IF REGENERATING THE DATA)

% Specify true values of the parameters
p1_true = 0.2;            % True fitness advantage for 1st vs 2nd
p2_true = 0.8;            % True fitness penalty for 1st vs 3rd
X0 = [0.2;0.8];           % Initial proportions

% Specify the data to which the parameters are being fit
N_pop = 1000;             % Population size for WF model
N_gen = 100;              % Number of generations in WF model (or t_end in replicator)
t_gen = 1;                % Time associated with each generation

% Specify parameter search ranges (to plot over)
p1_range = [0,1];
p2_range = [-0.5,1.5];
% Number of points to use in surface plotting (in each dimension)
Npts = 501;
% Number of points in between each plotted descent direction arrow
spacing = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATA REGENERATION

% Assume data is not to be regenerated unless necessary, or requested
if nargin < 1
    regenerate = false;
end

% Set data filename 
data_filename = 'DATA_FULL_fitSurfaces.mat';

% Define the type II replicator RHS for the two allele case
fitness_fun = @(x,p1,p2) [1 + p1, 1 + p1*p2; 1 + p1*p2, 1] * x;
df = @(x,p1,p2) x .* ( fitness_fun(x,p1,p2) ./ sum(x .* fitness_fun(x,p1,p2), 1) - 1 ) ;

% Generate the data if not present or asked to regenerate
if regenerate || ~exist(data_filename, 'file')
    
    % Prepare the timepoints/generations for observations
    tpts = 0:N_gen;
        
    % Data is a realisation from the Wright-Fisher model (type I selection)
    rng(1);
    X_data = wrightFisher(N_pop, N_gen, X0, @(x) fitness_fun(x,p1_true,p2_true), 2);
       
    % Approximate the derivatives using functional fit
    [F, Fdash] = constructSimplexFit(struct('t',tpts,'X',X_data));
    dX_data = Fdash(tpts);
        
    % Build the grid of points used to show the optimisation surface
    [p1m,p2m] = meshgrid( linspace(p1_range(1), p1_range(2), Npts), linspace(p2_range(1), p2_range(2), Npts) );

    % Evaluate the optimisation surfaces
    [Ny,Nx] = size(p1m);
    f_opt = zeros(Ny,Nx);
    df_opt = zeros(Ny,Nx);
    for i = 1:Ny
               
        if mod(i,10) == 0
            fprintf('Ran %g out of %g rows\n', i, Ny);
        end
        
        for j = 1:Nx
        
            % Run the replicator for these parameters and reshape
            [T,x_here] = runReplicator(2, t_gen, tpts, @(x) fitness_fun(x,p1m(i,j),p2m(i,j)), X0);
            x_here = x_here';
            % Evaluate RHS for derivatives
            dX_here = df(X_data, p1m(i,j), p2m(i,j) );
            % Calculate error
            f_opt(i,j) = norm(x_here - X_data,'fro');
            df_opt(i,j) = norm(dX_here - dX_data,'fro');
        
        end
    end

    % Save the data
    save(data_filename, 'p1_true', 'p2_true', 'p1m', 'p2m', 'f_opt', 'df_opt', 'X_data');

    % Safety pause
    pause(5);
    
end
    

%%% PLOTTING

% Specify the number of slots (controls spacing of figure panels)
N_slots = 11;

% Load in the data
load(data_filename, 'p1_true', 'p2_true', 'p1m', 'p2m', 'f_opt', 'df_opt');

% Load colormaps
load('extra_colormaps.mat','viridis');

% Prepare the figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
tiles_obj = tiledlayout(1,2*N_slots+1);

% First panel - nonlinear least squares
nexttile([1 N_slots]);
hold on;
surf(p1m,p2m,f_opt);
colormap(viridis);
shading flat;
view(2);
xlabel('Selection strength, $s$','Interpreter','LaTeX','Fontsize',24);
ylabel('Dominance, $h$','Interpreter','LaTeX','FontSize',24);
set(gca,'FontSize',24);
title('Nonlinear Least Squares','FontSize',28);
% Estimate gradients using the surface
gradx_f = diff(f_opt,1,2) ./ diff(p1m,1,2);
grady_f = diff(f_opt,1,1) ./ diff(p2m,1,1);
% Trim the extra row/column from each
gradx_f = gradx_f(1:end-1,:);
grady_f = grady_f(:,1:end-1);
% Convert to descent directions via normalisation
dec_xdir = -gradx_f ./ sqrt( gradx_f.^2 + grady_f.^2 );
dec_ydir = -grady_f ./ sqrt( gradx_f.^2 + grady_f.^2 );
% Prepare quiver plot data - plot only every Nth points
plotlocs = 1+round(spacing/2) : spacing : Npts - 1;
X = p1m(plotlocs,plotlocs) + p1m(1,2) - p1m(1,1);
Y = p2m(plotlocs,plotlocs) + p2m(2,1) - p2m(1,1);
Z = (max(f_opt(:)) + 1) * ones(size(X));   % 3D height is just "above surface" so it's seen
U = dec_xdir(plotlocs, plotlocs);
V = dec_ydir(plotlocs, plotlocs);
W = zeros(size(U));
% Add the quiver plot to the axes
quiver3(X,Y,Z,U,V,W, 0.5, 'LineWidth', 2, 'AutoScale','off','Color',[1 1 1],'ShowArrowHead','off','Marker','.','MarkerSize',25);
% Plot the true value
plot3(p1_true,p2_true,max(f_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.8, 0.0, 0.0],'MarkerFaceColor',[1.0, 0.2, 0.2]);
% Find and plot the minimum value
[~,minloc] = min(f_opt(:));
df_b1 = p1m(minloc);
df_c1 = p2m(minloc);
plot3(df_b1,df_c1,max(f_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.4, 0.4, 1.0],'MarkerFaceColor',[0.4, 1.0, 1.0]); 
% Ensure axes are correct
xlim([p1_range(1) p1_range(2)]);
ylim([p2_range(1) p2_range(2)]);
% Make axis square after all has been set
axis square;

% Gap panel
nexttile([1 1]);
set(gca,'Visible',false);

% Second panel - gradient regression
nexttile([1 N_slots]);
hold on;
surf(p1m,p2m,df_opt);
colormap(viridis);
shading flat;
view(2);
xlabel('Selection strength, $s$','Interpreter','LaTeX','Fontsize',24);
ylabel('Dominance, $h$','Interpreter','LaTeX','FontSize',24);
set(gca,'FontSize',24);
title('Gradient Matching','FontSize',28);
% Estimate gradients using the surface
gradx_f = diff(df_opt,1,2) ./ diff(p1m,1,2);
grady_f = diff(df_opt,1,1) ./ diff(p2m,1,1);
% Trim the extra row/column from each
gradx_f = gradx_f(1:end-1,:);
grady_f = grady_f(:,1:end-1);
% Convert to descent directions via normalisation
dec_xdir = -gradx_f ./ sqrt( gradx_f.^2 + grady_f.^2 );
dec_ydir = -grady_f ./ sqrt( gradx_f.^2 + grady_f.^2 );
% Prepare quiver plot data - plot only every Nth points
plotlocs = 1+round(spacing/2) : spacing : Npts - 1;
X = p1m(plotlocs,plotlocs) + p1m(1,2) - p1m(1,1);
Y = p2m(plotlocs,plotlocs) + p2m(2,1) - p2m(1,1);
Z = (max(df_opt(:)) + 1) * ones(size(X));   % 3D height is just "above surface" so it's seen
U = dec_xdir(plotlocs, plotlocs);
V = dec_ydir(plotlocs, plotlocs);
W = zeros(size(U));
% Add the quiver plot to the axes
quiver3(X,Y,Z,U,V,W, 0.5, 'LineWidth', 2, 'AutoScale','off','Color',[1 1 1],'ShowArrowHead','off','Marker','.','MarkerSize',25);
% Plot the true value
plot3(p1_true,p2_true,max(df_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.8, 0.0, 0.0],'MarkerFaceColor',[1.0, 0.2, 0.2]);
% Find and plot the minimum value
[~,minloc] = min(df_opt(:));
df_b1 = p1m(minloc);
df_c1 = p2m(minloc);
plot3(df_b1,df_c1,max(df_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.4, 0.4, 1.0],'MarkerFaceColor',[0.4, 1.0, 1.0]); 
% Ensure axes are correct
xlim([p1_range(1) p1_range(2)]);
ylim([p2_range(1) p2_range(2)]);
% Make axis square after all has been set
axis square;

% Tighten up the layout
tiles_obj.TileSpacing = 'compact';
tiles_obj.Padding = 'compact';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, X] = runReplicator(N_feat, t_gen, tpts, fitness, X0)

% Prepare a mass matrix to confine results to probability simplex
M = eye(N_feat);
M(N_feat,N_feat) = 0;

% Initialise ODE solve tolerance and most strict tolerance
tol = 1e-3;
min_tol = 1e-10;
running = true;
success = false;

% Disable warnings as we will try many tolerances
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

% Solve the replicator ODE
while running

    % Attempt solve using the current tolerance
    odeoptions = odeset('Mass',M,'AbsTol',tol,'RelTol',tol);
    [T,X] = ode15s( @(t,X) 1/t_gen * replicatorRHS(X,fitness,2,true), tpts, X0, odeoptions);

    % Terminate if integration successful or tolerance too low
    if abs(T(end) - tpts(end)) < 1e-10 || tol <= min_tol
        running = false;
        success = true;
    end

    % Terminate with failure if tolerance too small
    if tol <= min_tol && ~success
        running = false;
        fprintf('\n WARNING: There was a failure to solve the ODE even with strictest tolerance! Output will be truncated.\n');
    end

    % Otherwise, decrease tolerance and try again
    tol = tol * 0.01;

end

% Re-enable warnings to not bother the user's session
warning('on','MATLAB:ode15s:IntegrationTolNotMet');
        