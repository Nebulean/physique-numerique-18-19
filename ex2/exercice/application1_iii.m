clear

% on fait la simulation
cmd = './Exercice2 configuration.in nsteps=500 E=6e4 B0=3 vx0=0 vy0=4e5 tfin=1.09321e-7 schema=RK2 output=app1_iii.out';
disp(cmd);
system(cmd);

vE = 2e4;

% on traite les données
d = load('app1_iii.out');
t = d(:,1);
x = d(:,2);
y = d(:,3);
energy = d(:,6);

% On applique la transformation galiléenne sur la position
x = x - vE.*t;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

% On plot
fig1 = figure;
p = plot(x,y);

set(p, 'LineWidth',1.5);
set(gca, 'fontsize',20);

xlabel('x [m]');
ylabel('y [m]');

grid on;

fig2 = figure;
p = plot(t, energy, 'LineWidth',1.5);
set(gca,'fontsize',20);
xlabel('t [s]');
ylabel('E [J]');

grid on;


saveas(fig1, 'graphs/app1_iii_traj','epsc');
saveas(fig2, 'graphs/app1_iii_ene','epsc');