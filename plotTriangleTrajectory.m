function plotTriangleTrajectory(X_WF, X_rep, figtitle, show_legend)
% This function uses the triangle plotting function to plot the provided
% Wright-Fisher trajectory and associated replicator trajectory on the
% probability simplex (a triangle)

% Create figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
hold on;

% Plot Wright-Fisher data
trianglePlot(NaN, X_WF', 'k.', 'MarkerSize', 25);

% Plot replicator data
trianglePlot(NaN, X_rep', 'r', 'LineWidth',3);

% Create legend if requested
if show_legend
   
    % Specify number of things to be plotted
    % (for easy modification of code later if needed)
    N_data = 2;
    % Prepare legend text
    leg_txt = cell(1,N_data*7);
    legend_names = {'Wright-Fisher', 'Replicator'};
    for k = 1:2
        
        % Set location in legend text array
        loc = (k-1)*7;
        % Add six blank entries (3 edges + 3 vertices of the triangle)
        for m = 1:6
            leg_txt{loc+m} = '';
        end
        % Now add the specified label
        leg_txt{loc+7} = legend_names{k};
        
    end
    
    % Display legend
    legend(leg_txt,'Location','NorthWest','FontSize',20,'box','off');
    
end

% Append title
title(figtitle,'FontSize',24,'Interpreter','LaTeX');

end
