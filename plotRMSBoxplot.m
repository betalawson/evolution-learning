function plotRMSBoxplot(simResults, x_txts, title_txt, fig_filename, save_figs, plot_pop)
%
% This function plots the RMSE data from the input structure (stored in the
% fields 'l2_mat', as well as 'l2_perfect' if 'plot_pop' is provided a
% value of true). The user provides the x tick labels as 'x_txts',
% additional title text as 'title_txt', and also the figure's filename as
% 'fig_filename' that will be used if the variable 'save_figs' is set to a
% value of true.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the colours to use in plotting the two selection types
plot_clrs = [   [ 0.6, 0.0, 0.0 ];          % Dark red
                [ 0.0, 0.0, 0.6 ];          % Dark blue
                [ 1.0, 0.4, 0.4 ];          % Light red
                [ 0.4, 0.4, 1.0 ]  ];       % Light blue
          
% Specify the fraction of data to show on the box plot whiskers
keep_frac = 0.95;

% Figure dimensions
fig_x = 0.3;
fig_dx = 0.25;
fig_y = 0.3;
fig_dy = fig_dx * (1920/1080);
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume perfect value will be plotted
if nargin < 6
    plot_pop = false;
end

% Take logarithms of the data because this is an error plot
l2_mat = log10( simResults.l2_mat );
l2_perfect = log10( simResults.l2_perfect );



%%% PLOTTING
 
% Prepare the figure
fig_obj = figure('units','normalized','Position',[fig_x fig_y fig_dx fig_dy]);
hold on;
      
% Trim down the data to only the requested confidence interval
l2_mat = trimToCI(l2_mat, keep_frac);
    
% Plot the boxplot showing the main data
boxplot_obj = boxplot2(l2_mat, 'Whisker', Inf);
    
% Add colours to the plot
for ii = 1:size(l2_mat,2)
    structfun(@(x) set(x(ii,:), 'color',plot_clrs(ii,:), 'markeredgecolor',plot_clrs(ii,:), 'LineWidth',2), boxplot_obj);
end
    
% Add to the plot the parameter learned for perfect data (if requested)
% Axis limits are also set here
if plot_pop
    for ii = 1:size(l2_mat,2)
        plot(Nx-0.25+0.1*ii,l2_perfect(ii),'.', 'MarkerSize',40, 'MarkerEdgeColor',plot_clrs(ii,:));
    end
    axis([0.5 Nx+0.5  -0.1+min([l2_mat(:);l2_perfect(:)]) 0.1+max([l2_mat(:); l2_perfect(:)]) ]);
else
    axis([0.5 Nx+0.5  -0.1+min(l2_mat(:)) 0.1+max(l2_mat(:)) ]);
end
        
% Other plot clean-up
xticks(1:Nx);
xticklabels(x_txts);
set(gca, 'FontSize',20, 'LineWidth',2);
if plot_pop
    xlabel('$N$', 'FontSize',24, 'Interpreter','latex');
end
ylabel('$\log_{10} \mathrm{RMSE}$', 'FontSize',24, 'Interpreter','latex');
title(['Trajectory $L_2$ Error (',title_txt,')'],'FontSize',20,'Interpreter','latex');
    
% Save the figure and then close it (if save flag set)
if save_figs
    saveas(fig_obj, [fig_filename,'.eps'], 'epsc');
    close(fig_obj);
end

end