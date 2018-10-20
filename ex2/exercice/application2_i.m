% Nom du fichier d'output a analyser
f = {};
f{1} = 'App2pos.out';
f{2} = 'App2neg.out';

fig = {};

% output = {};
% output{1} = 'graphs/graphApp2pos';
% output{2} = 'graphs/graphApp2neg';

% proton

filename = f{1};
    
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

clear output;

% Figures

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

fig1 = figure;
p1 = plot(x,y);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');

% fig2 = figure;
% p2 = plot(vx,vy);
% set(gca, 'fontsize',16);
% set(p1, 'LineWidth',1.5);
% axis equal;
% grid on;
% xlabel('$v_x [m/s]$');
% ylabel('$v_y [m/s]$');

fig2 = figure;
p1 = plot(t,mu);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('$t [s]$');
ylabel('$\mu [J/T]$');
legend('$B(x)$')

fig3 = figure;
p1 = plot(t,x,t,y);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('x,y [m]');
legend('x','y');

fig4 = figure;
p1 = plot(t,vx,t,vy);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('$v_x,v_y [m/s]$');
legend('$v_x$','Proton in magnetic field $B(x)$');

fig5 = figure;
p1 = plot(t,energy);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('E [J]');
legend('Proton');

% antiproton

filename = f{2};
    
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

clear output;

% Figures

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

fig6 = figure;
p1 = plot(x,y);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');

% fig7 = figure;
% p2 = plot(vx,vy);
% set(gca, 'fontsize',16);
% set(p1, 'LineWidth',1.5);
% axis equal;
% grid on;
% xlabel('$v_x [m/s]$');
% ylabel('$v_y [m/s]$');

fig8 = figure;
p1 = plot(t,x,t,y);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('x,y [m]');
legend('x','y');

fig9 = figure;
p1 = plot(t,vx,t,vy);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('$v_x,v_y [m/s]$');
legend('$v_x$','Antiproton in magnetic field $B(x)$');

fig10 = figure;
p1 = plot(t,energy);
set(gca, 'fontsize',20);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('E [J]');
legend('Antiproton');

saveas(fig1, 'graphs/app2_ipos_traj', 'epsc');
saveas(fig2, 'graphs/app2_ipos_mu', 'epsc');
saveas(fig3, 'graphs/app2_ipos_xyt', 'epsc');
saveas(fig4, 'graphs/app2_ipos_vxvyt', 'epsc');
saveas(fig5, 'graphs/app2_ipos_ene', 'epsc');

saveas(fig6, 'graphs/app2_ineg_traj', 'epsc');
saveas(fig7, 'graphs/app2_ineg_vit', 'epsc');
saveas(fig8, 'graphs/app2_ineg_xyt', 'epsc');
saveas(fig9, 'graphs/app2_ineg_vxvyt', 'epsc');
saveas(fig10, 'graphs/app2_ineg_ene', 'epsc');

