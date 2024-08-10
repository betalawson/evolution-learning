function trianglePlot(t, v, varargin)
% This function takes as input a series of vectors in the probability
% simplex for R3, and plots them as a trajectory through this simplex as
% visualised by a triangle where each vertex corresponds to a pure
% population of one type, (1,0,0), (0,1,0) or (0,0,1). This is naturally
% achieved using the barycentric co-ordinates for the triangle.

% Set up the vertices of an equilateral triangle with side length one, and
% centroid at the origin:
r1 = [-1/2,-sqrt(3)/6];
r2 = [1/2,-sqrt(3)/6];
r3 = [0,2*sqrt(3)/6];

% Transform the input set of vectors into positions on this triangle
rv = v(:,1) * r1 + v(:,2) * r2 + v(:,3) * r3;

% Plot vertices as black dots
plot(r1(1),r1(2),'k.','MarkerSize',40);
plot(r2(1),r2(2),'k.','MarkerSize',40);
plot(r3(1),r3(2),'k.','MarkerSize',40);
% Plot black lines between these
plot([r1(1),r2(1)],[r1(2),r2(2)],'k','LineWidth',3);
plot([r1(1),r3(1)],[r1(2),r3(2)],'k','LineWidth',3);
plot([r2(1),r3(1)],[r2(2),r3(2)],'k','LineWidth',3);
% Plot the path on the triangle
plot(rv(:,1),rv(:,2),varargin{:});

% Adjust plot axis limits (adds a fixed margin in longest direction and square axis)
margin = 0.02;
axis equal;
xlim([r1(1)-margin, r2(1)+margin]);
ylim([r1(2)-margin, r3(2)+margin]);

% Turn off the ticks and hide axis
yticks([]);
xticks([]);
ax = gca;
ax.XAxis.Visible = false;
ax.YAxis.Visible = false;

end

