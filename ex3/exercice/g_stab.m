%% Définition des constantes
g=9.81;
d=0.035;
L=0.01;
Omega=sqrt(2*g)/d;
theta0=pi-0.01;

%% On lance la simulation
cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=%f kappa=0.1 theta0=%0.15f thetadot0=0 dt=0.02 tFin=20 output=g_stab.out", Omega, d, theta0)
system(cmd);
disp(cmd);

%% ANALYSE
g=load("g_stab.out");

t=g(:,1);
theta=g(:,2);

f1=figure;
plot(t, theta);
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);
grid on;
xlabel("$t[s]$");
ylabel("$\theta [rad]$");

saveas(f1, "graphs/g_stabk", 'epsc');