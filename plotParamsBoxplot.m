function plotParamsBoxplot(simResults, x_txts, xlabel_txt, title_txt, fig_filename, save_figs, add_legend, plot_pop)
%
% This function plots the values of the inferred parameters from the input
% structure (stored in the fields 'params_mat' and 'params_true', as well
% as 'params_perfect' if 'plot_pop' is provided a value of true). The user
% provides the x tick labels as 'x_txts', additional title text as
% 'title_txt', and also the figure's filename as 'fig_filename' that will
% be used if the variable 'save_figs' is set to a value of true.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the colours to use in plotting the two selection types
plot_clrs = [   [ 0.6, 0.0, 0.0 ];          % Dark red
                [ 1.0, 0.4, 0.4 ];          % Light red                
                [ 0.0, 0.0, 0.6 ];          % Dark blue
                [ 0.4, 0.4, 1.0 ]  ];       % Light blue
          
% Specify how these should be labelled in (optional) legend
legend_txts = {'LS S_I', 'GM S_I', 'LS S_{II}', 'GM S_{II}'};
            
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

% Assume legend text will not be added unless specifically requested
if nargin < 7
    add_legend = false;
end

% Assume perfect value will not be plotted if not specifically requested
if nargin < 8
    plot_pop = false;
end

% Extract the parameter names and count how many
param_names = fieldnames( simResults.params_mat{1,1,1} );
N_params = length(param_names);

% Count the number of x ticks provided
Nx = length(x_txts);


%%% PLOTTING FOR EACH PARAMETER

% Loop over each parameter
for p = 1:N_params

    % Grab out this parameter's name
    Pname = param_names{p};
    
    % Extract the 3D array of values for the current parameter
    p_mat = cellfun(@(x) x.(Pname), simResults.params_mat);
        
    % Extract the true value for the parameter
    p_true = simResults.params_true.(Pname);
    
    % Prepare the figure
    fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]);
    hold on;
    
    % Plot the true value first
    plot([-1, Nx+1], [p_true p_true], 'k', 'LineWidth', 1.5);
    
    % Trim down the data to only the requested confidence interval
    p_mat = trimToCI(p_mat, keep_frac);
    
    % Re-order the data to put type I and type II together
    p_mat = p_mat(:,[1,3,2,4],:);
    
    % Plot the boxplot showing the main data
    boxplot_obj = boxplot2(p_mat, 'Whisker', Inf);
    
    % Add colours to the plot
    for ii = 1:size(p_mat,2)
        structfun(@(x) set(x(ii,:), 'color',plot_clrs(ii,:), 'markeredgecolor',plot_clrs(ii,:), 'LineWidth',2), boxplot_obj);
    end
    
    % Add to the plot the parameter learned for perfect data (if requested)
    if plot_pop
        
        % Extract the perfect values for the current parameter
        p_perfect = cellfun(@(x) x.(Pname), simResults.params_perfect); 
        % Re-order the perfect data to put type I and type II together
        p_perfect = p_perfect([1,3,2,4]);
        % Plot these in their appropriate colours
        for ii = 1:size(p_mat,2)
            plot(Nx-0.25+0.1*ii,p_perfect(ii),'.', 'MarkerSize',40, 'MarkerEdgeColor',plot_clrs(ii,:));
        end
        
    end
    
    % Set up the plot axis limits
    y_absmin = p_true - y_abslength/2;
    y_absmax = p_true + y_abslength/2;
    pmin = min(p_mat(:));
    pmax = max(p_mat(:));
    ymin = max([y_absmin; pmin-0.1*(pmax-pmin)]);
    ymax = min([y_absmax; pmax+0.1*(pmax-pmin)]);
    axis([0.5 Nx+0.5 ymin ymax]);
    
    % Other plot clean-up
    xticks(1:Nx);
    xticklabels(x_txts);
    set(gca, 'FontSize',20, 'LineWidth',2);
    xlabel(xlabel_txt, 'FontSize',24, 'Interpreter','latex');
    ylabel(['$',Pname,'$'], 'FontSize',24, 'Interpreter','latex');
    if ~isempty(title_txt)
        title(['Estimated $',Pname,'$ (',title_txt,')'], 'FontSize',20, 'Interpreter','latex');
    else
        title(['Estimated $',Pname,'$'], 'FontSize',20, 'Interpreter','latex');
    end
    
    % Add a legend if requested
    if add_legend
        for ii = 1:size(p_mat,2)
            text(Nx-0.5, ymax - (ymax - ymin)*(0.05 + 0.075*(ii-1)), legend_txts{ii}, 'FontSize', 16, 'FontName', 'DejaVu Sans Mono', 'Color', plot_clrs(ii,:));
        end
    end
    
    % Save the figure and then close it (if save flag set)
    if save_figs
        saveas(fig_obj, [fig_filename,'_',Pname,'.eps'], 'epsc');
        close(fig_obj);
    end

end