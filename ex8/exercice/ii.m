%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
delta = 0;
x0 = 0;
sigma_norm = 0.06;
n = 14;
tfin = 5000;
Ninters = 300;
m = 1.;
h = 1.;

dt = 1;

k0=2*pi*n/(xR-xL);

cmd = sprintf("./Exercice8 configuration.in output=ii xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", xL, xR, omega, delta, x0, sigma_norm, n, tfin, Ninters, dt);


%% Simulations
disp(cmd);
system(cmd);

%% Data
data=load('ii_obs.out');
t=data(:,1);
xmoy = data(:,5);
pmoy = data(:,7);

% xth=x0*cos(omega.*t)+sqrt(2*h/(m*omega))*sin(omega.*t);
% pth=m*(-x0*sin(omega.*t)+sqrt(2*h*omega/m)*cos(omega.*t));

xth=x0*cos(omega.*t)+k0/omega*sin(omega.*t);
pth=-x0*m*omega*sin(omega.*t)+k0*cos(omega.*t);

xdiff=abs(xmoy-xth);
pdiff=abs(pmoy-pth);

[xmax,xid]=findpeaks(xdiff);
[pmax,pid]=findpeaks(pdiff);

%% Figures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

% f1=figure;
% hold on;
% set(gca, 'fontsize', 25);
% set(gca, 'LineWidth',1.5);
% plot(t,xmoy, 'linewidth',1.2);
% plot(t,xth, 'linewidth',1.2);
% xlabel('$t ~ [t_P]$');
% ylabel('$x ~ [\ell_P]$');
% grid on;
% box on;
% legend('$\langle x \rangle (t)$', 'Newtonian $x(t)$');
% hold off;
% 
% f2=figure;
% hold on;
% set(gca, 'fontsize', 25);
% set(gca, 'LineWidth',1.5);
% plot(t,pmoy,'linewidth',1.2);
% plot(t,pth,'linewidth',1.2);
% xlabel('$t ~ [t_P]$');
% ylabel('$p ~ [m_P ~ c]$');
% grid on;
% box on;
% legend('$\langle p \rangle (t)$', 'Newtonian $p(t)$');
% hold off;
% 
f3=figure;
hold on;
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);
plot(t,xmoy-xth, 'linewidth',1.2);
xlabel('$t ~ [t_P]$');
ylabel('$\langle x \rangle (t)-x(t)~ [\ell_P]$');
grid on;
box on;
hold off;

f4=figure;
hold on;
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);
plot(t,pmoy-pth,'linewidth',1.2);
xlabel('$t ~ [t_P]$');
ylabel('$\langle p \rangle (t)-p(t)~ [m_P ~ c]$');
grid on;
box on;
hold off;

f5=figure;
hold on;
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);
plot(t(xid),xmax, 'x', 'markersize',10, 'linewidth',1.2);
px=polyfit(t(xid),xmax,1);
plot(t(xid),px(1)*t(xid)+px(2), 'linewidth',1.2);
xlabel('$t ~ [t_P]$');
ylabel('$|\langle x \rangle (t)-x(t)| ~ [\ell_P]$');
grid on;
box on;
legend('maxima of $|\langle x \rangle (t)-x(t)|$', "fit: slope="+num2str(px(1)),'location','northwest');
hold off;

f6=figure;
hold on;
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);
plot(t(pid),pmax,'x', 'markersize',10, 'linewidth',1.2);
pp=polyfit(t(pid),pmax,1);
plot(t(pid),pp(1)*t(pid)+pp(2), 'linewidth',1.2);
xlabel('$t ~ [t_P]$');
ylabel('$|\langle p \rangle (t)-p(t)|~ [m_P ~ c]$');
grid on;
box on;
legend('maxima of $|\langle p \rangle (t)-p(t)|$', "fit: slope="+num2str(pp(1)),'location','northwest');
hold off;

% saveas(f1,'graphs/ii_x','epsc');
% saveas(f2,'graphs/ii_p','epsc');
saveas(f3,'graphs/ii_diffx','epsc');
saveas(f4,'graphs/ii_diffp','epsc');
saveas(f5,'graphs/ii_diffxsl','epsc');
saveas(f6,'graphs/ii_diffpsl','epsc');