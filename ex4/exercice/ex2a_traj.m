%% On génère les données
dtad="true";
epsilon=1e-5
atm="true";
tFin=172800;
cmd = sprintf("./Exercice4 configuration.in tFin=%f atm=%s dt=8 dtad=%s epsilon=%f output=ex2a_traj.out", tFin, atm, dtad, epsilon);
system(cmd)
%% On load les données
d = load("ex2a_traj.out");

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

ax3 = d(:,15);
ay3 = d(:,16);

P = d(:,17);

%% Trajectoire
fig1=figure
hold on;
% on plot les endroits initiaux
pbaspect([1 1 1]);
daspect([1 1 1]);
plot(x3(1), y3(1), 'x', 'Color','green');
centerOfEarth = [x1(1), y1(1)];
plotCircle(centerOfEarth, 6371000, 500, 'blue');

% puis on plot les positions.
plot(x1(2:end), y1(2:end), '-', 'Color', 'blue', 'LineWidth', 1.2);
plot(x3(2:end), y3(2:end), '-', 'Color', 'green', 'LineWidth', 1.2);

grid on;
hold off;

% nsteps = length(t)

%% Etude de convergence de l'acceleration
% fig2=figure

%% Quelques calculs supplémentaires
maxAccel = max( sqrt(ax3.^2 + ay3.^2) );
g = 9.80665; % source: https://fr.wikipedia.org/wiki/G_(acc%C3%A9l%C3%A9ration)
fprintf("L'acceleration max de Apollo est %0.2f [m/s^2], soit %0.2fG.\n", maxAccel, maxAccel/g);

maxPower = max( abs( P ) );
fprintf("La puissance max de la force de frottement est %0.2f [W].\n", maxPower);


function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb);

    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
       
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end