function plotRMSEBoxplot(simResults, x_txts, xlabel_txt, title_txt, fig_filename, save_figs, add_legend, plot_pop)
%
% This function plots the RMSE data from the input structure (stored in the
% fields 'l2_mat', as well as 'l2_perfect' if 'plot_pop' is provided a
% value of true). The user provides the x tick labels as 'x_txts',
% additional title text as 'title_txt', and also the figure's filename as
% 'fig_filename' that will be used if the variable 'save_figs' is set to a
% value of true.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the colours to use in plotting the two selection types and two
% regression methods
plot_clrs = [   [ 0.6, 0.0, 0.0 ];          % Dark red
                [ 1.0, 0.4, 0.4 ];          % Light red
                [ 0.0, 0.0, 0.6 ];          % Dark blue
                [ 0.4, 0.4, 1.0 ]  ];       % Light blue
            
% Specify how these should be labelled in (optional) legend
legend_txts = {'LS S_I', 'GM S_I', 'LS S_{II}', 'GM S_{II}'};
          
% Specify the fraction of data to show on the box plot whiskers
keep_frac = 1;

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

% Count the number of x labels given
Nx = length(x_txts);


%%% PLOTTING
 
% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]);
hold on;
      
% Take logarithms of the data because this is an error plot
l2_mat = log10( simResults.l2_mat );

% Trim down the data to only the requested confidence interval
l2_mat = trimToCI(l2_mat, keep_frac);   

% Re-order the columns to put type I results and type II results together
l2_mat = l2_mat(:,[1,3,2,4],:);

% Plot the boxplot showing the main data
boxplot_obj = boxplot2(l2_mat, 'Whisker', Inf);

% Add colours to the plot
for ii = 1:size(l2_mat,2)
    structfun(@(x) set(x(ii,:), 'color',plot_clrs(ii,:), 'markeredgecolor',plot_clrs(ii,:), 'LineWidth',2), boxplot_obj);
end
    
% Add to the plot the parameter learned for perfect data (if requested)
% Axis limits are also set here
if plot_pop
    
    % Read out data for the perfect case
    l2_perfect = log10( simResults.l2_perfect );
    
    % Similarly re-order this data
    l2_perfect = l2_perfect([1,3,2,4]);
    
    % Add to the plot
    for ii = 1:size(l2_mat,2)
        plot(Nx-0.25+0.1*ii,l2_perfect(ii),'.', 'MarkerSize',40, 'MarkerEdgeColor',plot_clrs(ii,:));
    end
    
    % Grab out minimum and maximum values for the y-axis
    ymin = min([l2_mat(:); l2_perfect(:)]) - 0.1;
    ymax = max([l2_mat(:); l2_perfect(:)]) + 0.1;
        
    % Specify the axes to fit all content on the plot
    axis([0.5 Nx+0.5 ymin ymax]);
    
else
    
    % Grab out minimum and maximum values for the y-axis
    ymin = min(l2_mat(:)) - 0.1;
    ymax = max(l2_mat(:)) + 0.1;
    
    % If not plotting the perfect result, specify axes without extra space
    axis([0.5 Nx+0.5 ymin ymax]);
    
end
        
% Other plot clean-up
xticks(1:Nx);
xticklabels(x_txts);
yticks( floor(ymin*2)/2:0.5:ceil(ymax*2)/2 );
set(gca, 'FontSize',20, 'LineWidth',2);
xlabel(xlabel_txt, 'FontSize',24, 'Interpreter','latex');
ylabel('$\log_{10} \mathrm{RMSE}$', 'FontSize',24, 'Interpreter','latex');
if ~isempty(title_txt)
    title(title_txt,'FontSize',20,'Interpreter','latex');
else
    title('Trajectory $L_2$ Error','FontSize',20,'Interpreter','latex');
end

% Add a legend if requested
if add_legend
    for ii = 1:size(l2_mat,2)
        text(Nx-0.5, 0.1 + max(l2_mat(:)) - (0.2 + max(l2_mat(:)) - min(l2_mat(:)))*(0.05 + 0.075*(ii-1)), legend_txts{ii}, 'FontSize', 16, 'FontName', 'DejaVu Sans Mono', 'Color', plot_clrs(ii,:));
    end 
end
    
% Save the figure and then close it (if save flag set)
if save_figs
    saveas(fig_obj, [fig_filename,'.eps'], 'epsc');
    close(fig_obj);
end

end