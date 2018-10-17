% Nom du fichier d'output a analyser
f = {};
f{1} = 'App2pos.out';
f{2} = 'App2neg.out';

fig = {};

% output = {};
% output{1} = 'graphs/graphApp2pos';
% output{2} = 'graphs/graphApp2neg';

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
set(gca, 'fontsize',16);
set(p1, 'LineWidth',1.5);
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');

fig2 = figure;
p2 = plot(vx,vy);
set(gca, 'fontsize',16);
set(p1, 'LineWidth',1.5);
axis equal;
grid on;
xlabel('$v_x [m/s]$');
ylabel('$v_y [m/s]$');

fig3 = figure;
plot(t,x,t,y);
set(gca, 'fontsize',16);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('x,y [m]');
legend('x','y');

fig4 = figure;
plot(t,vx,t,vy);
set(gca, 'fontsize',16);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('$v_x,v_y [m/s]$');
legend('$v_x$','$v_y$');

fig5 = figure;
plot(t,energy);
set(gca, 'fontsize',16);
set(p1, 'LineWidth',1.5);
grid on;
xlabel('t [s]');
ylabel('E [J]');

saveas(fig1, 'graphs/app2_i_traj', 'epsc');
saveas(fig2, 'graphs/app2_i_vit', 'epsc');
saveas(fig3, 'graphs/app2_i_xyt', 'epsc');
saveas(fig4, 'graphs/app2_i_vxvyt', 'epsc');
saveas(fig5, 'graphs/app2_i_ene', 'epsc');

