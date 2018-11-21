g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);

Omega = omega0;

% They're watching.
% kappa=0.0005; d=0.04; theta0=2*pi/3; thetadot0=0.;

% 1984
% kappa=0.; d=0.04; theta0=2*pi/3; thetadot0=0.;

% Bernd das Brot
% kappa=0.0005; d=0.04; theta0=2*pi/3; thetadot0=0.; Omega=3*Omega

% Solar System. (Homoculus ?)
% kappa=0.0001; d=0.001; theta0=2*pi/3; thetadot0=10; Omega=Omega;

% rainbow
kappa=0.0001; d=0.001; theta0=pi-1e-8; thetadot0=10; Omega=1.7*Omega;

n = 100; % NE PAS TOUCHER
tfin = 30000*2*pi/Omega; % NE PAS TOUCHER
dt = 2*pi/(n*Omega);

% couleurs
purple = [0.5, 0, 1];
blue = [0, 0, 1];
cyan = [0, 1, 1];
green = [0, 1, 0];
yellow = [1, 1, 0];
orange = [1, 0.5, 0];
red = [1, 0, 0];

color = {purple, blue, cyan, green, yellow, orange, red};
legend = {}

fig2=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

legend = [];

for i=1:7
    Omega=(1.45 + (i-1)*0.035)*omega0;
    
    cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f, d=%0.15f kappa=%0.15f theta0=%0.15f thetadot0=%0.15f tFin=%0.15f dt=0.02 sampling=%d output=g_funfunfun_rainbow%d.out", Omega, d, kappa, theta0, thetadot0, tfin, n, i)
    
    
    
    legend{i} = sprintf("%0.2f", Omega);
    %system(cmd);
    
    toLoad = sprintf("g_funfunfun_rainbow%d.out", i);
    data = load(toLoad);
    theta = wrapToPi(data(:,2));
    thetadot = data(:,3);
    plot(theta, thetadot, '.', 'Color', color{i});
end

grid on;

legend(legend);

xlabel("$\theta$ [rad]");
ylabel("$\dot{\theta}$ [rad/s]");

hold off


% cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f, d=%0.15f kappa=%0.15f theta0=%0.15f thetadot0=%0.15f tFin=%0.15f dt=0.02 sampling=%d output=g_funfunfun.out", Omega, d, kappa, theta0, thetadot0, tfin, n);
% system(cmd)

%% On load les données
d = load("g_funfunfun.out");

theta = d(:,2);
thetadot = d(:,3);

theta=wrapToPi(theta);
%% On plot
fig=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(theta, thetadot, '.');

xlabel("$\theta$ [rad]");
ylabel("$\dot{\theta}$ [rad/s]");

grid on;

hold off;