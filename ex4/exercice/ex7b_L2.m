%% On load les données
cmd = "./Exercice4 configLagrange.in x3=4.446515204924e8 y3=0 vx3=0 vy3=-1183.48732826579"

system(cmd);
d = load("Lagrange.out");

Omega = 2.66160639e-6;

t = d(:,1);
dt = d(:,2);

x1 = d(:,3);
y1 = d(:,4);

x2 = d(:,7);
y2 = d(:,8);

x3 = d(:,11);
y3 = d(:,12);

vx3 = d(:,13);
vy3 = d(:,14);

dAT = sqrt((x3-x1).^2 + (y3-y1).^2);
dAL = sqrt((x3-x2).^2 + (y3-y2).^2);

old = x1;
x1 = x1.*cos(Omega*t) - y1.*sin(Omega*t);
y1 = old.*sin(Omega*t) + y1.*cos(Omega*t);

old = x2;
x2 = x2.*cos(Omega*t) - y2.*sin(Omega*t);
y2 = old.*sin(Omega*t) + y2.*cos(Omega*t);

old = x3;
x3 = x3.*cos(Omega*t) - y3.*sin(Omega*t);
y3 = old.*sin(Omega*t) + y3.*cos(Omega*t);

old = vx3;
vx3 = vx3.*cos(Omega*t) - vy3.*sin(Omega*t);
vy3 = old.*sin(Omega*t) + vy3.*cos(Omega*t);

%% On plot les endroits initiaux.
fig=figure
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);


pbaspect([1 1 1]);
daspect([1 1 1]);
plot(x3(1), y3(1), 'x', 'Color','magenta','MarkerSize',10);
text(x3(1) - 0.3e8, 0.4e8, "L2", 'Color', 'magenta');
centerOfEarth = [x1(end), y1(end)];
plotCircle(centerOfEarth, 6371000, 500, 'blue');
centerOfMoon = [x2(end), y2(end)];
plotCircle(centerOfMoon, 1737500, 50, 'red');

% puis on plot les positions.
plot(x1, y1, '-', 'Color', 'blue', 'LineWidth', 1.2);
plot(x2, y2, '-', 'Color', 'red', 'LineWidth', 1.2);

v=sqrt(vx3(2:end).^2 + vy3(2:end).^2);
coloredLinePlot(x3(2:end), y3(2:end), v, "Speed $[m/s]$");

xlabel("x [m]");
ylabel("y [m]");

grid on;
box on;

hold off;

% saveas(fig, 'graphs/ex7b_L1_traj','epsc');
print(fig, "graphs/ex7b_L2_traj", '-dpng', '-r500');


f2=figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(t,dAT, t,dAL);

xlabel("Time t [s]");
ylabel("Distance [m]")
grid on;
box on;
legend('Distance to Earth','Distance to Moon');


hold off;

saveas(f2, 'graphs/ex7b_L2_dist','epsc');


nsteps = length(t)

function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb);

    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
       
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end


function c=coloredLinePlot(x, y, col, label)
    %Source: https://blogs.mathworks.com/videos/2014/08/12/coloring-a-line-based-on-height-gradient-or-some-other-value-in-matlab/
    surface('XData', [x x],             ... % N.B.  XYZC Data must have at least 2 cols
        'YData', [y y],             ...
        'ZData', zeros(numel(x),2), ...
        'CData', [col col],             ...
        'FaceColor', 'none',        ...
        'EdgeColor', 'interp',      ...
        'Marker', 'none');
    
    c=colorbar();
    colormap jet;	
    c.Label.Interpreter = 'latex';
    c.Label.String = label;
end