function plotAllParamsBoxplot(simResults, pop_levels, selection_type, title_txt, fig_filename, save_figs, add_legend, legend_loc)
%
% This function plots the values of the inferred parameters from the input
% structure (stored in the fields 'params_mat' and 'params_true', all on a
% single plot, for two provided population levels.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the colours to use in plotting the two selection types
plot_clrs = [   [ 0.6, 0.0, 0.0 ];          % Dark red                
                [ 0.0, 0.0, 0.6 ];          % Dark blue
                [ 1.0, 0.4, 0.4 ];          % Light red
                [ 0.4, 0.4, 1.0 ]  ];       % Light blue
            
% Specify the fraction of data to show on the box plot whiskers
keep_frac = 1;

% Specify the total maximum length of the y-axis (in parameter space)
y_abslength = 2.5;

% Figure dimensions
fig_x = 0.3;
fig_dx = 0.25;
fig_y = 0.3;
fig_dy = fig_dx * (1920/1080);
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INITIALISATION

% Assume legend text will not be added unless specifically requested
if nargin < 7
    add_legend = false;
end
% Assume legend will appear in top right unless specifically requested
if nargin < 8
    legend_loc = 'top_right';
end

% Count out how many population levels are requested
Np = length(pop_levels);

% Prepare the legend texts (if one will be added)
if add_legend
        
    % Legend text is method name, plus population size
    for n = 1:Np
        legend_txts{2*n-1} = ['LS ', num2str( simResults.pop_sizes(pop_levels(n)) )];
        legend_txts{2*n} = ['GM ', num2str( simResults.pop_sizes(pop_levels(n)) )];
    end
    
end


%%% DATA PREPARATION

% Extract the parameter names and count how many
param_names = fieldnames( simResults.params_mat{1,1,1} );
N_params = length(param_names);

% Initialise the new data structure
all_mat = NaN( N_params, 2*Np, size(simResults.params_mat,3) );

% Loop over each parameter, and fill out the new data structure
for p = 1:N_params
    
    % Grab out this parameter's name
    Pname = param_names{p};
    
    % Extract the 3D array of values for the current parameter
    p_mat = cellfun(@(x) x.(Pname), simResults.params_mat);
    
    % Grab out the data for each requested population level
    for n = 1:Np
        % Least squares has type I/II in slots [1,2]
        all_mat(p,2*n-1,:) = p_mat(pop_levels(n),selection_type,:);
        % Gradient matching has type I/II in slots [3,4]
        all_mat(p,2*n,:) = p_mat(pop_levels(n),2+selection_type,:);
    end
    
    % Convert all data for this parameter to be the difference as compared
    % to the true value
    all_mat(p,:,:) = all_mat(p,:,:) - simResults.params_true.(Pname);
    
    %%% Parameter name cleanup
    
    % Find the first character(s)
    char_find = regexp(Pname,'[a-zA-Z]+','match');
    char_txt = char_find{1};
    % Find the numeric character(s)
    num_find = regexp(Pname,'[0-9]+','match');
    num_txt = num_find{1};
    % Convert any selection for heterozygotes to 'H' - notation in paper
    if length(num_txt) == 2 && strcmp(char_txt,'s')
        char_txt = 'H';
    end
    % Now set up the LaTeX math version of parameter for use on the plot:
    %    $ <name>_{ <numbers> } $
    plot_paramnames{p} = ['$',char_txt,'_{',num_txt,'}$'];
    
end


%%% PLOTTING

% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]);
hold on;

% Plot the true value first - zero because we subtracted it already
plot([-1, N_params+1], [0 0], 'k', 'LineWidth', 1.5);
    
% Trim down the data to only the requested confidence interval
all_mat = trimToCI(all_mat, keep_frac);
    
% Plot the boxplot showing the main data
boxplot_obj = boxplot2(all_mat, 'Whisker', Inf);
    
% Add colours to the plot
for ii = 1:size(p_mat,2)
    structfun(@(x) set(x(ii,:), 'color',plot_clrs(ii,:), 'markeredgecolor',plot_clrs(ii,:), 'LineWidth',2), boxplot_obj);
end
    
    
% Set up the plot axis limits
y_absmin = - y_abslength/2;
y_absmax = + y_abslength/2;
pmin = min(all_mat(:));
pmax = max(all_mat(:));
ymin = max([y_absmin; pmin-0.1*(pmax-pmin)]);
ymax = min([y_absmax; pmax+0.1*(pmax-pmin)]);
axis([0.5 N_params+0.5 ymin ymax]);

% Other plot clean-up
xticks(1:N_params);
xticklabels(plot_paramnames);
set(gca,'TickLabelInterpreter','LaTeX');
set(gca, 'FontSize',20, 'LineWidth',2);
ylabel('Error in Estimate', 'FontSize',24, 'Interpreter','latex');

% Prepare type text for title
if selection_type == 1
    type_txt = 'Type I';
else
    type_txt = 'Type II';
end
% Add title
if ~isempty(title_txt)
    title([title_txt, ' (',type_txt,')'], 'FontSize',20, 'Interpreter','latex');
else
    title(['Parameter Error (',type_txt,')'], 'FontSize',20, 'Interpreter','latex');
end

% Add a legend if requested
if add_legend
    switch lower(legend_loc)
        case 'top_right'
            for ii = 1:length(legend_txts)
                text(N_params-1, ymax - (ymax - ymin)*(0.05 + 0.075*(ii-1)), legend_txts{ii}, 'FontSize', 16, 'FontName', 'DejaVu Sans Mono', 'Color', plot_clrs(ii,:));
            end
        case 'top_left'
            for ii = 1:length(legend_txts)
                text(0.75, ymax - (ymax - ymin)*(0.05 + 0.075*(ii-1)), legend_txts{ii}, 'FontSize', 16, 'FontName', 'DejaVu Sans Mono', 'Color', plot_clrs(ii,:));
            end
        case 'bottom_right'
            for ii = 1:length(legend_txts)
                text(N_params-1.35, ymin + (ymax - ymin) * (0.05 + 0.075*(length(legend_txts)-1) - 0.075*(ii-1) ), legend_txts{ii}, 'FontSize', 16, 'FontName', 'DejaVu Sans Mono', 'Color', plot_clrs(ii,:));
            end
        case 'bottom_left'
            for ii = 1:length(legend_txts)
                text(0.75, ymin + (ymax - ymin) * (0.05 + 0.075*(length(legend_txts)-1) - 0.075*(ii-1) ), legend_txts{ii}, 'FontSize', 16, 'FontName', 'DejaVu Sans Mono', 'Color', plot_clrs(ii,:));
            end
    end
end

% Save the figure and then close it (if save flag set)
if save_figs
    saveas(fig_obj, [fig_filename,'.eps'], 'epsc');
    close(fig_obj);
end

end