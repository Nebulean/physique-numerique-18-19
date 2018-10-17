% on fait la simulation
cmd = './Exercice2 configuration.in nsteps=500 E=6e4 B0=3 vx0=-2e4 vy0=4e5 tfin=1.09321e-7 schema=EulerCromer output=app1_iii.out';
disp(cmd);
system(cmd);

% on traite les données
d = load('app1_iii.out');
x = d(:,2);
y = d(:,3);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

% On plot
fig = figure;
p = plot(x,y);

set(p, 'LineWidth',1.5);
set(gca, 'fontsize',16);

xlabel('x [m]');
ylabel('y [m]');