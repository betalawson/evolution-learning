function FIGURE2_lossLandscapes(regenerate)
% This function plots the loss surface to be optimised in fitting the
% parameters of a basic bi-allelic inheritance problem to synthetic data.
% This is compared with the corresponding loss surface for gradient
% matching the corresponding derivative data.


%%% PARAMETER SETTINGS (ONLY USED IF REGENERATING THE DATA)

% Specify true values of the parameters
s_true = 0.2;             % Selective advantage
h_true = 0.8;             % Dominance in heterozygotes
X0 = [0.2;0.8];           % Initial proportions

% Specify the data to which the parameters are being fit
data_WF = false;           % Use Wright-Fisher model or replicator
N_pop = 1000;             % Population size for WF model
N_gen = 50;               % Number of generations in WF model (or t_end in replicator)
N_obs = 11;               % Number of observation points in the data
t_gen = 1;                % Time associated with each generation

% Specify parameter search ranges (to plot over)
s_range = [0,1];
h_range = [-0.5, 1.5];
% Number of points to use in surface plotting (in each dimension)
Npts = 501;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATA REGENERATION

% Assume data is not to be regenerated unless necessary, or requested
if nargin < 1
    regenerate = false;
end

% Set data filename according to whether using perfect replicator data or
% generated Wright-Fisher data
if data_WF
    data_filename = 'DATA_WF_lossSurfaces.mat';
else
    data_filename = 'DATA_REP_lossSurfaces.mat';
end

% Generate the data if not present or asked to regenerate
if regenerate || ~exist(data_filename, 'file')
   
    % Define the type I replicator RHS for the two allele case
    fitness = @(x,s,h) [1 + s, 1 + h*s; 1 + h*s, 1] * x;
    df = @(x,s,h) x .* ( fitness(x,s,h) - sum(x .* fitness(x,s,h), 1) );
    
    % Prepare the timepoints/generations for observations
    tpts = round(linspace(0,N_gen,N_obs));
    
    % Generate synthetic data including derivative
    if data_WF
        
        % Set seed for repeatability
        rng(7);
        
        % Simulate the Wright-Fisher model with type I selection
        X_data = wrightFisher(N_pop, N_gen, X0, @(x) fitness(x,s_true,h_true), 1);
        
        % Cut down to only the observable locations
        X_data = X_data(:,tpts+1);
        
        % Approximate the derivatives using functional fit
        [F, Fdash] = constructSimplexFit(struct('t',tpts,'X',X_data));
        dX_data = Fdash(tpts);
        
    else

        % Run the replicator with the true fitness
        [T,X_data] = runReplicator(2, t_gen, tpts, @(x) fitness(x,s_true,h_true), X0);
        % Put array in the right shape
        X_data = X_data';
        % Generate associated derivative data direct from replicator
        dX_data = df(X_data, s_true, h_true);
        
    end

    % Build the grid of points used to show the optimisation surface
    [sm,hm] = meshgrid( linspace(s_range(1), s_range(2), Npts), linspace(h_range(1), h_range(2), Npts) );

    % Evaluate the optimisation surfaces
    [Ny,Nx] = size(sm);
    f_opt = zeros(Ny,Nx);
    df_opt = zeros(Ny,Nx);
    for i = 1:Ny
        parfor j = 1:Nx
        
            % Run the replicator for these parameters and reshape
            [T,x_here] = runReplicator(2, t_gen, tpts, @(x) fitness(x,sm(i,j),hm(i,j)), X0);
            x_here = x_here';
            % Evaluate RHS for derivatives
            dX_here = df(X_data, sm(i,j), hm(i,j) );
            % Calculate error
            f_opt(i,j) = norm(x_here - X_data,'fro');
            df_opt(i,j) = norm(dX_here - dX_data,'fro');
        
        end
    end

    % Save the data
    save(data_filename, 's_true', 'h_true', 'sm', 'hm', 'f_opt', 'df_opt');

    % Safety pause
    pause(5);
    
end
    

%%% PLOTTING

% Specify the particle start locations for the gradient descent plots
p0 = [0.1, -0.4;
      0.9, -0.4;
      0.1, 1.4;
      0.9, 1.4];

% Specify the number of slots (controls spacing of figure panels)
N_slots = 11;

% Load in the data
load(data_filename, 's_true', 'h_true', 'sm', 'hm', 'f_opt', 'df_opt');

% Load colormaps
load('extra_colormaps.mat','viridis');

% Prepare the figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
tiles_obj = tiledlayout(1,2*N_slots+1);

% Left panel - standard regression
nexttile([1 N_slots]);
hold on;
surf(sm,hm,f_opt);
colormap(viridis);
shading flat;
view(2);
xlabel('Selection Advantage, s');
ylabel('Dominance, h');
set(gca,'FontSize',24);
title('Nonlinear Least Squares','FontSize',28);
% Add trajectories to the plot using gradient descent function, below
addTrajectories(p0,sm,hm,f_opt,s_range,h_range,s_true,h_true);
% Plot the true value        
plot3(s_true,h_true,max(f_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.8, 0.0, 0.0],'MarkerFaceColor',[1.0, 0.2, 0.2]);

% "Middle" panel - gap
nexttile([1 1]);
set(gca,'Visible',false);

% Right panel - gradient regression
nexttile([1 N_slots]);
hold on;
surf(sm,hm,df_opt);
colormap(viridis);
shading flat;
view(2);
xlabel('Selection Advantage, s');
ylabel('Dominance, h');
set(gca,'FontSize',24);
title('Gradient Matching','FontSize',28);
% Add trajectories to the plot using gradient descent function, below
addTrajectories(p0,sm,hm,df_opt,s_range,h_range,s_true,h_true);
% Plot the true value
plot3(s_true,h_true,max(df_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.8, 0.0, 0.0],'MarkerFaceColor',[1.0, 0.2, 0.2]);

% Tighten up the layout
tiles_obj.TileSpacing = 'compact';
tiles_obj.Padding = 'compact';



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
    [T,X] = ode15s( @(t,X) 1/t_gen * replicatorRHS(X,fitness,1,true), tpts, X0, odeoptions);

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



function addTrajectories(p0,sm,hm,f,s_range,h_range,s_true,h_true)

% Specify learning rate for gradient descent
learn_rate = 0.01;

% Specify tolerance for gradient descent
tol = 1e-4;

% Specify maximum number of gradient descent steps
max_iters = 100000;

% Create gradient descent trajectories
for k = 1:size(p0,1)
    
    % Initialise the particle at this location
    part_s = p0(k,1);
    part_h = p0(k,2);
    part_traj = [part_s; part_h];
    
    % Loop until sufficiently close to the minimum
    running = true; iters = 0;
    while running
        
        % Read out the four corners surrounding the particle's current location
        sloc = find( sm(1,:) > part_s, 1);
        hloc = find( hm(:,1) > part_h, 1);
        % Corner locations
        s_p = sm(1,sloc);
        s_m = sm(1,sloc-1);
        h_p = hm(hloc,1);
        h_m = hm(hloc-1,1);
        % Get function values in these corners
        f_pp = f(hloc,sloc);
        f_mp = f(hloc,sloc-1);
        f_pm = f(hloc-1,sloc);
        f_mm = f(hloc-1,sloc-1);
        
        % Get the gradient using bilinear interpolation
        grad_s = ( (part_h - h_m)*(f_pp - f_mp) + (h_p - part_h)*(f_pm - f_mm) ) / ( (s_p - s_m) * (h_p - h_m) );
        grad_h = ( (part_s - s_m)*(f_pp - f_pm) + (s_p - part_s)*(f_mp - f_mm) ) / ( (s_p - s_m) * (h_p - h_m) );
        
        % Take a gradient descent step and add it to trajectory
        part_s = part_s - learn_rate * grad_s;
        part_h = part_h - learn_rate * grad_h;
        part_traj(:,end+1) = [part_s; part_h];
        
        % Check for closeness
        if sqrt( (part_s - s_true)^2 + (part_h - h_true)^2 ) < tol
            running = false;
        end
        
        % Check for leaving the plot area
        if part_s < s_range(1) || part_s > s_range(2) || part_h < h_range(1) || part_h > h_range(2)
            running = false;
        end
        
        % Check for iteration count
        iters = iters + 1;
        if iters >= max_iters
            running = false;
        end
        
    end
    
    %plot3(part_traj(1,:),part_traj(2,:),max(f(:))+1 * ones(1,size(part_traj,2)),'Color',[1 1 1],'LineWidth',4);
    
    % Add the starting point as a dot
    %plot3(part_traj(1,1),part_traj(2,1),max(f(:))+1,'.','Color',[1 1 1],'MarkerSize',35);
    plot3(part_traj(1,1),part_traj(2,1),max(f(:))+1,'.','Color',[1 1 1],'MarkerSize',20);
    
    % Add this trajectory to the plot with two colours of different thickness
    
    plot3(part_traj(1,:),part_traj(2,:),max(f(:))+1 * ones(1,size(part_traj,2)),'Color',[1 1 1],'LineWidth',2.5);
    

    
end