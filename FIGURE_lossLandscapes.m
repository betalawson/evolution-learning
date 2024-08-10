function FIGURE_lossLandscapes(regenerate)
% This function plots the loss surface to be optimised in fitting the
% parameters of a basic bi-allelic inheritance problem to synthetic data.
% This is compared with the corresponding loss surface for gradient
% matching the corresponding derivative data.


%%% PARAMETER SETTINGS (ONLY USED IF REGENERATING THE DATA)

% Specify true values of the parameters
s_true = 0.5;             % Selective advantage
h_true = 0.8;             % Dominance in heterozygotes
x0 = [0.2;0.8];           % Initial proportions

% Specify parameter search ranges (to plot over)
s_range = [0,2];
h_range = [-1, 2];
% Number of points to use in surface plotting (in each dimension)
Npts = 501;

% Specify the timepoints to observe
Nobs = 26;
t_end = 50;
tpts = linspace(0,t_end,Nobs);


%%% DATA REGENERATION

% Assume data is not to be regenerated unless necessary, or requested
if nargin < 1
    regenerate = false;
end

% Generate the data if not present or asked to regenerate
if regenerate || ~exist('DATA_lossSurfaces.mat', 'file')
   
    % Define the replicator RHS for two allele case
    fitness = @(x,s,h) [1 + s, 1 + h*s; 1 + h*s, 1] * x;
    df = @(x,s,h) x .* ( fitness(x,s,h) - sum(x .* fitness(x,s,h), 1) );

    % Generate synthetic data including derivative
    odeoptions = odeset('AbsTol',1e-8,'RelTol',1e-8);
    [T,x_data] = ode15s( @(t,x) df(x,s_true,h_true), tpts, x0, odeoptions );
    x_data = x_data';
    dx_data = df(x_data, s_true, h_true);

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
            dx_here = df(x_data, sm(i,j), hm(i,j) );
            % Calculate error
            f_opt(i,j) = norm(x_here - x_data,'fro');
            df_opt(i,j) = norm(dx_here - dx_data,'fro');
        
        end
    end

    % Save the data
    save('DATA_lossSurfaces.mat', 's_true', 'h_true', 'sm', 'hm', 'f_opt', 'df_opt');

    % Safety pause
    pause(5);
    
end
    

%%% PLOTTING

% Load in the data
load('DATA_lossSurfaces.mat', 's_true', 'h_true', 'sm', 'hm', 'f_opt', 'df_opt');

% Plot the data
figure; subplot(1,2,1); surf(sm,hm,f_opt); view(2); shading flat; subplot(1,2,2); surf(sm,hm,df_opt); view(2); shading flat;