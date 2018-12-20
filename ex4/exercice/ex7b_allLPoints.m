% L1
x1 = 3.219948720257e8;
y1 = 0;

% L2
x2 = 4.446515204924e8;
y2 = 0;

% L3
x3 = -3.866963974114e8;
y3 = 0;

% L4
x4 = 187697755.46;
y4 = 333201542.1;

% L5
x5 = 187697755.46;
y5 = -333201542.1;


% Earth
ex = -4.676244537e6;
ey = 0;

% Moon
mx = 380071755.5;
my = 0;

fig=figure
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);


pbaspect([1 1 1]);
daspect([1 1 1]);

xlim([-4.8e8, 4.8e8]);
ylim([-4.8e8, 4.8e8]);

color = {'blue', 'green', 'red', 'cyan', 'magenta'};

plot(x1, y1, 'x', 'MarkerSize', 10, 'Color', color{1});
plot(x2, y2, 'x', 'MarkerSize', 10, 'Color', color{2});
plot(x3, y3, 'x', 'MarkerSize', 10, 'Color', color{3});
plot(x4, y4, 'x', 'MarkerSize', 10, 'Color', color{4});
plot(x5, y5, 'x', 'MarkerSize', 10, 'Color', color{5});

text(x1 - 0.15e8, y1 + 0.3e8, "L1", 'Color', color{1});
text(x2 - 0.15e8, y2 + 0.3e8, "L2", 'Color', color{2});
text(x3 + 0.25e8, y3 + 0.3e8, "L3", 'Color', color{3});
text(x4 - 0.15e8, y4 + 0.35e8, "L4", 'Color', color{4});
text(x5 - 0.15e8, y5 - 0.35e8, "L5", 'Color', color{5});

%centerOfEarth = [ex, ey];
%plotCircle(centerOfEarth, 6371000, 500, 'blue');
plot(ex, ey, 'o', 'MarkerSize', 10);

plot(mx, my, 'o', 'MarkerSize', 10);

% Toutes les lignes
%L2 - L3
x = [x2, x3];
y = [y2, y3];
fit = poly_approx(x, y, 1, 2);
plot(fit(1,:), fit(2,:), '-', 'Color','black');

%m2 - L4
x = [mx x4];
y = [my y4];
fit = poly_approx(x, y, 1, 2);
plot(fit(1,:), fit(2,:), '-', 'Color','black');

%m2 - L5
x = [mx x5];
y = [my y5];
fit = poly_approx(x, y, 1, 2);
plot(fit(1,:), fit(2,:), '-', 'Color','black');

%m1 - L4
x = [ex x4];
y = [ey y4];
fit = poly_approx(x, y, 1, 2);
plot(fit(1,:), fit(2,:), '-', 'Color','black');

%m1 - L5
x = [ex x5];
y = [ey y5];
fit = poly_approx(x, y, 1, 2);
plot(fit(1,:), fit(2,:), '-', 'Color','black');

% Le cercle autour
% radius = abs(ex + mx);
radius = sqrt((ex - x4)^2 + (ey - y4)^2);
% radius = abs(ex + mx)/2;
centerOfEarth = [ex, 0];
plotCircle(centerOfEarth, radius, 500, 'black');


%lgd = ["L1", "L2", "L3", "L4", "L5", "Earth", "Moon"]
%legend(lgd);

xlabel("x [m]");
ylabel("y [m]");

grid on;
box on;

hold off;


saveas(fig, 'graphs/ex6b_allLPoints','epsc');



function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb)

    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
       
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end

function [polynome, slope] = poly_approx(x, y, ordre, steps)
    pf = polyfit(x, y, ordre);
    slope = pf(1);
    T = linspace(min(x), max(x), steps);

    n = ordre + 1;

    polynome = zeros(2,length(T));
    for i=1:n
       polynome(2,:) = polynome(2,:) + pf(i)*T.^(n-i);
    end

    polynome(1,:) = T;
end