% Nom du fichier d'output a analyser
filename = 'app1_ii.out';

    
% Chargement des donnees
output = load(filename);

% Extraction des quantites d'interet
t = output(:,1);
x = output(:,2);
y = output(:,3);
vx = output(:,4);
vy = output(:,5);
energy = output(:,6);
mu = output(:,7);

clear output

% Figures

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

fig1 = figure;
p1 = plot(x,y);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('x [m]');
ylabel('y [m]');

fig2 = figure;
p2 = plot(vx,vy);
set(gca, 'fontsize',20);
set(p2, 'LineWidth',1.5);
grid on;
xlabel('$v_x$ [m/s]');
ylabel('$v_y$ [m/s]');

fig3 = figure;
p3 = plot(t,energy);
set(gca, 'fontsize',20);
set(p3, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('E [J]');
%set(gca, 'XScale', 'log');
%set(gca,'YScale', 'log');

saveas(fig1, 'graphs/app1_ii_traj', 'epsc');
saveas(fig2, 'graphs/app1_ii_vit', 'epsc');
saveas(fig3, 'graphs/app1_ii_ene', 'epsc');