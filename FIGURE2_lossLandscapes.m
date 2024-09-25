function FIGURE2_lossLandscapes(regenerate)
% This function plots the loss surface to be optimised in fitting the
% parameters of a basic bi-allelic inheritance problem to synthetic data.
% This is compared with the corresponding loss surface for gradient
% matching the corresponding derivative data.


%%% PARAMETER SETTINGS (ONLY USED IF REGENERATING THE DATA)

% Specify true values of the parameters
s_true = 0.2;             % Selective advantage
h_true = 0.8;             % Dominance in heterozygotes
x0 = [0.2;0.8];           % Initial proportions

% Specify the data to which the parameters are being fit
data_WF = true;           % Use Wright-Fisher model or replicator
N_pop = 5000;             % Population size for WF model
N_gen = 50;               % Number of generations in WF model (or t_end in replicator)
N_obs = 11;               % Number of observation points in the data

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
if regenerate || ~exist(data_filename', 'file')
   

    
    % Define the type I replicator RHS for the two allele case
    fitness = @(x,s,h) [1 + s, 1 + h*s; 1 + h*s, 1] * x;
    df = @(x,s,h) x .* ( fitness(x,s,h) - sum(x .* fitness(x,s,h), 1) );
    % Prepare ODE solution options
    odeoptions = odeset('AbsTol',1e-8,'RelTol',1e-8);
    
    % Prepare the timepoints/generations for observations
    tpts = round(linspace(0,N_gen,N_obs));
    
    % Generate synthetic data including derivative
    if data_WF
        
        % Set seed for repeatability
        rng(7);
        
        % Simulate the Wright-Fisher model with type I selection
        X_data = wrightFisher(N_pop, N_gen, x0, fitness, 1);
        
        % Cut down to only the observable locations
        X_data = X_data(:,tpts);
        
        % Approximate the derivatives using functional fit
        [F, Fdash] = constructSimplexFit(struct('t',tpts,'X',X_data));
        dX_data = Fdash(tpts);
        
    else
        
        % Simulate the replicator
        [T,X_data] = ode15s( @(t,x) df(x,s_true,h_true), tpts, x0, odeoptions );
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
        for j = 1:Nx
        
            % Solve ODE for these parameters
            [T,x_here] = ode15s( @(t,x) df(x, sm(i,j), hm(i,j)), tpts, x0, odeoptions );
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
hold on;
surf(sm,hm,f_opt);
colormap(viridis);
plot3(s_true,h_true,max(f_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[1.0, 0.2, 0.2],'MarkerFaceColor',[1.0, 0.2, 0.2]);
shading flat;
view(2);
xlabel('Selection Advantage, s');
ylabel('Dominance, h');
set(gca,'FontSize',24);
title('Standard Regression','FontSize',28);

% "Middle" panel - gap
nexttile([1 1]);
set(gca,'Visible',false);

% Right panel - gradient regression
nexttile([1 N_slots]);
hold on;
surf(sm,hm,df_opt);
colormap(viridis);
plot3(s_true,h_true,max(df_opt(:))+1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[1.0, 0.2, 0.2],'MarkerFaceColor',[1.0, 0.2, 0.2]);
shading flat;
view(2);
xlabel('Selection Advantage, s');
ylabel('Dominance, h');
set(gca,'FontSize',24);
title('Gradient Matching Regression','FontSize',28);

% Tighten up the layout
tiles_obj.TileSpacing = 'compact';
tiles_obj.Padding = 'compact';
