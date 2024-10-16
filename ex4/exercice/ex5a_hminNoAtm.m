%% Initialisation

dtad="true";
epsilon=1e-4;
atm="false";
tFin=172800;

% On va appliquer une rotation légère sur la vitesse initiale d'Apollo
nsimul = 25;
center=-.309676
range = 1*10^(-6)
angle = linspace(center-range, center+range, nsimul);

% On importe la vitesse initiale
vx       = -1178.49845516659;
vy       = -597.5766552;

% On calcul les nouvelles vitesses
newvx = zeros(nsimul, 1);
newvy = zeros(nsimul, 1);

newvx = vx*cos(angle) - vy*sin(angle);
newvy = vx*sin(angle) + vy*cos(angle);

%% SIMULATIONS
for i=1:nsimul
    cmd = sprintf("./Exercice4 conf3corps.in tFin=%f atm=%s dt=8 dtad=%s epsilon=%f output=ex5a_hminNoAtm%d.out vx3=%0.15f vy3=%0.15f", tFin, atm, dtad, epsilon, i, newvx(i), newvy(i));
    disp(cmd)
    system(cmd)
end

%% On charge et traite ces données et on plot la trajectoire de chaque simulation.
f1=figure
hold on;
pbaspect([1 1 1]);
daspect([1 1 1]);

% hmin
RT = 6378.1 * 1000; % earth's radius
distTH = 10e3 + RT;
mindist = zeros(nsimul, 1);

for i=1:nsimul
    d = load(sprintf("ex5a_hminNoAtm%d.out",i));
    
    t = d(:,1);
    dt = d(:,2);
    
    x1 = d(:,3);
    y1 = d(:,4);

    x2 = d(:,7);
    y2 = d(:,8);

    x3 = d(:,11);
    y3 = d(:,12);
    
    %% On calcul la distance minimale
    % Distance numérique avec les dt
    dist = sqrt((x1 - x3).^2 + (y1 - y3).^2);
    % Premièrement, on calcule le point où il y a la distance minimale.
    [tmp, index] = min(dist);
    if index+1 > length(dist)
        fit = polyfit(t(index-2:index), dist(index-2:index), 2);
    elseif index-1 < 1
        fit = polyfit(t(index:index+2), dist(index:index+2), 2);
    else
        fit = polyfit(t(index-1:index+1), dist(index-1:index+1), 2);
    end
    A = fit(1); B = fit(2); C = fit(3);
    
    % On calcule le minimum à l'aide de l'analyse
    mindist(i, 1) = abs(C - B^2/(4*A) - distTH); %TODO: Retirer les abs lorsque le résultat sera correct

    % puis on plot les positions.
    plot(x1(2:end), y1(2:end), '-', 'Color', 'blue');
    % plot(x2(2:end), y2(2:end), '-');
    plot(x3(2:end), y3(2:end), '-');

    grid on;
    
end

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

% on plot les endroits initiaux
plot(x3(1), y3(1), 'x', 'Color','green');
centerOfEarth = [x1(end), y1(end)];
plotCircle(centerOfEarth, RT, 500, 'blue');
set(gca, 'fontsize', 22);
xlabel("x [m]");
ylabel("y [m]");
box on;
hold off;

% on plot hmin
f2=figure;
plot (angle,mindist,'x');
set(gca, 'fontsize', 22);
xlabel('Angle difference with initial velocity direction [rad]');
ylabel('Difference of minimal height w/ $h_{min}$ [m]');
box on;
grid on;

saveas(f1,'graphs/ex5a_traj.eps','epsc');
saveas(f2,'graphs/ex5a_hmin.eps','epsc');

function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb);

    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
       
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end