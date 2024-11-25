function FIGURE9b_realDatalossLandscapes(regenerate)
% This function plots the loss surface to be optimised in fitting the
% parameters of a basic bi-allelic inheritance problem to synthetic data.
% This is compared with the corresponding loss surface for gradient
% matching the corresponding derivative data.


%%% PARAMETER SETTINGS (ONLY USED IF REGENERATING THE DATA)

% Load in the panaxia data
load('panaxia_data.mat','panaxia_data');
% Time values are years, shifted to start at zero
data.t = panaxia_data.Year - panaxia_data.Year(1);
% X values are the frequency of the unfavoured allele, then favoured
data.X = [panaxia_data.MedionigraAlleleFreq'; panaxia_data.RegularAlleleFreq'];

% Set up the problem object associated
t_gen = 1;
X0 = data.X(:,1);

% Specify parameter search ranges (to plot over)
s_range = [-1,0.5];
h_range = [-0.25,1.25];
% Number of points to use in surface plotting (in each dimension)
Npts = 251;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATA REGENERATION

% Assume data is not to be regenerated unless necessary, or requested
if nargin < 1
    regenerate = false;
end

% Set data filename according to whether using perfect replicator data or
% generated Wright-Fisher data
data_filename = 'DATA_REAL_lossSurfaces.mat';


% Generate the data if not present or asked to regenerate
if regenerate || ~exist(data_filename, 'file')
   
    % Define the type II replicator RHS for the two allele case
    fitness = @(x,s,h) [1 + s, 1 + h*s; 1 + h*s, 1] * x;
    df = @(x,s,h) x .* ( fitness(x,s,h) ./ sum(x .* fitness(x,s,h), 1) - 1);
    
    % Prepare the timepoints/generations for observations
    tpts = data.t;
    
    % Approximate the derivatives using functional fit
    [F, Fdash] = fitSimplexConstrainedFunction(data);
    dX_data = Fdash(tpts);
    
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
            dX_here = df(data.X, sm(i,j), hm(i,j) );
            % Calculate error, including shrinkage
            f_opt(i,j) = norm(x_here - data.X,'fro');
            df_opt(i,j) = norm(dX_here - dX_data,'fro');
                        
            %clf;
            
            %subplot(1,2,1); hold on;
            %plot(data.t,x_here(1,:),'r.','MarkerSize',25);
            %plot(data.t,data.X(1,:),'k.','MarkerSize',25);
            %title(['Error is',num2str(norm(x_here - data.X,'fro'))]);
            
            %subplot(1,2,2); hold on;
            %plot(data.X(1,:),dX_here(1,:),'r.','MarkerSize',25);
            %plot(data.X(1,:),dX_data(1,:),'k.','MarkerSize',25);
            %title(['Error is',num2str(norm(dX_here - dX_data,'fro'))]);
            %pause(0.1);
            
        end
    end

    % Save the data
    save(data_filename, 'sm', 'hm', 'f_opt', 'df_opt');

    % Safety pause
    pause(5);
    
end
   
%%% EQUATION LEARNING

% Prepare default options
ELoptions = addDefaultOptions(struct('symmetric_payoff',true),2);

% Find parameters using least squares
library = evolutionLearning(data, 2, setfield(ELoptions, 'regression_type', 'ls'));
A = libraryToPayoff(library,2);
LSparams = extractBialleleParams(A,2);

% Find parameters using least squares
[library, fitfun] = evolutionLearning(data, 2, setfield(ELoptions, 'regression_type', 'grad'));
A = libraryToPayoff(library,2);
ELparams = extractBialleleParams(A,2);


%%% PLOTTING

% Specify the number of slots (controls spacing of figure panels)
N_slots = 101;

% Load in the data
load(data_filename, 'sm', 'hm', 'f_opt', 'df_opt');

% Prepare the figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
tiles_obj = tiledlayout(1,2*N_slots+10);

% Left panel - standard regression
nexttile([1 N_slots]);
hold on;
surf(sm,hm,log10(f_opt/min(f_opt(:))));
colormap(flipud(bone));
caxis([0,0.25]);
shading flat;
view(2);
xlabel('Selection Advantage, s');
ylabel('Dominance, h');
set(gca,'FontSize',24);
title('Nonlinear Least Squares','FontSize',28);
axis square;

% Plot the hs = constant guideline
s_line = linspace(s_range(1),s_range(2),1001);
h_line = LSparams.hs ./ s_line;
h_line(h_line > 1+h_range(2)) = NaN;
h_line(h_line < -1+h_range(1)) = NaN;
plot3(s_line, h_line, ones(length(s_line),1), 'r', 'LineWidth', 2.5);
% Plot the value found using standard regression
plot3(LSparams.s, LSparams.hs / LSparams.s, 1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.8, 0.0, 0.0],'MarkerFaceColor',[1.0, 0.2, 0.2]);


% Limit axis
xlim([s_range(1), s_range(2)]);
ylim([h_range(1), h_range(2)]);

% "Middle" panel - gap
nexttile([1 1]);
set(gca,'Visible',false);

% Right panel - gradient regression
nexttile([1 N_slots]);
hold on;
surf(sm,hm,log10(df_opt/min(df_opt(:))));
colormap(flipud(bone));
caxis([0,0.25]);
shading flat;
view(2);
xlabel('Selection Advantage, s');
ylabel('Dominance, h');
set(gca,'FontSize',24);
title('Gradient Matching','FontSize',28);
axis square;

% Plot the hs = constant guideline
s_line = linspace(s_range(1),s_range(2),101);
h_line = ELparams.hs ./ s_line;
h_line(h_line > 1+h_range(2)) = NaN;
h_line(h_line < -1+h_range(1)) = NaN;
plot3(s_line, h_line, ones(length(s_line),1), 'r', 'LineWidth', 2.5);
% Plot the value found using gradient matching regression
plot3(ELparams.s, ELparams.hs / ELparams.s, 1,'pentagram','MarkerSize',25,'MarkerEdgeColor',[0.8, 0.0, 0.0],'MarkerFaceColor',[1.0, 0.2, 0.2]);

% Limit axis
xlim([s_range(1), s_range(2)]);
ylim([h_range(1), h_range(2)]);

% Limit axis
xlim([s_range(1), s_range(2)]);
ylim([h_range(1), h_range(2)]);

% Add colourbar
clrbar_obj = colorbar();
clrbar_obj.Position = [0.94, 0.145, 0.01, 0.75];
title(clrbar_obj,'log_{10}(L/L_{min})');
clr_ticklabels = clrbar_obj.TickLabels;
clr_ticklabels{end} = ['\geq',clr_ticklabels{end}];
clrbar_obj.TickLabels = clr_ticklabels;


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