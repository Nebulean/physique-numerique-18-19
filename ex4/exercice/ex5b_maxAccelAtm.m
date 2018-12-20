%% Initialisation

dtad="true";
epsilon=1e-4;
atm="true";
tFin=172800;

% On va appliquer une rotation légère sur la vitesse initiale d'Apollo
nsimul = 20;
center=-0.309177
range = 1*10^(-5)
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
    cmd = sprintf("./Exercice4 conf3corps.in tFin=%f atm=%s dt=8 dtad=%s epsilon=%f output=ex5a_maxAccelAtm%d.out vx3=%0.15f vy3=%0.15f", tFin, atm, dtad, epsilon, i, newvx(i), newvy(i));
    disp(cmd)
    system(cmd)
end

%% On charge et traite ces données et on plot la trajectoire de chaque simulation.
fig_traj=figure
hold on;
pbaspect([1 1 1]);
daspect([1 1 1]);
maxAccel = zeros(nsimul, 1);
RT = 6378.1 * 1000; % earth's radius
g = 9.80665;

for i=1:nsimul
    d = load(sprintf("ex5a_maxAccelAtm%d.out",i));
    
    t = d(:,1);
    dt = d(:,2);
    
    x1 = d(:,3);
    y1 = d(:,4);

    x2 = d(:,7);
    y2 = d(:,8);

    x3 = d(:,11);
    y3 = d(:,12);
    
    ax3 = d(:,15);
    ay3 = d(:,16);
    
    maxAccel(i,1) = max( sqrt(ax3.^2 + ay3.^2) ) / g;


    % puis on plot les positions.
    plot(x1(2:end), y1(2:end), '-', 'Color', 'blue');
    % plot(x2(2:end), y2(2:end), '-');
    plot(x3(2:end), y3(2:end), '-');

    grid on;
    
end

% on plot les endroits initiaux
plot(x3(1), y3(1), 'x', 'Color','green');
centerOfEarth = [x1(end), y1(end)];
plotCircle(centerOfEarth, RT, 500, 'blue');
hold off;

% on plot l'accel
fig = figure;
hold on;
plot(angle, maxAccel, 'x');
grid on;
hold off;

function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb);

    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
       
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end