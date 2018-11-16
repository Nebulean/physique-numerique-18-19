clear
%% Définition des constantes
g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);

Omega = 2*omega0;

n = 500; % On choisit le nombre de de pas de temps que l'on veut avoir dans chaque période.

tfin = 200*2*pi/Omega; % Le cas où i = 10000 est celui de base, mais j'ai pris un peu plus.

dt = 2*pi/(n*Omega);

%% CHOIX DE DIFFERENTES CONDITIONS INITIALES

% CHAOS, CHAOS
theta01 = 5*pi/6;
thetadot0 = 0.;
figname = "chaos";

% NOT CHAOS, BORING
% theta01 = 1e-6;
% thetadot0 = 1.;
% figname = "notchaos";

% deuxième graph
theta02 = theta01+1e-8;

%% SIMULATIOOOOONS

cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.05 kappa=0.1 theta0=%0.15f thetadot0=%0.15f dt=%0.15f tFin=%f sampling=%d output=f_first.out", Omega, theta01, thetadot0, dt, tfin, n)
system(cmd);
disp(cmd);

cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.05 kappa=0.1 theta0=%0.15f thetadot0=%0.15f dt=%0.15f tFin=%f sampling=%d output=f_second.out", Omega, theta02, thetadot0, dt, tfin, n)
system(cmd);
disp(cmd);

%% ANALYSE

d1 = load("f_first.out");
t = d1(:,1);
theta1 = d1(:,2);
thetadot1 = d1(:,3);

d2 = load("f_second.out");
theta2 = d2(:,2);
thetadot2 = d2(:,3);

%différence
d = sqrt((thetadot1-thetadot2).^2 + Omega*Omega*(theta1-theta2).^2);

%% PLOTILDE
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

f1=figure;
hold on;
theta1=wrapToPi(theta1);
theta2=wrapToPi(theta2);
plot(theta1, thetadot1, '.', theta2, thetadot2, '.');
set(gca, 'fontsize', 22);
xlabel("$\theta$ [rad]");
ylabel("$\dot{\theta}$ [rad/s]");
legend("$\theta_0=\pi$","$\theta_0=\pi+1 \cdot 10^{-8}$");
grid on;
hold off;

f2=figure;
hold on;
set(gca, 'YScale', 'log');
plot(t,d);
set(gca, 'fontsize', 22);
grid on
hold off;

wheretosave = sprintf("graphs/f_sens_%s", figname);
saveas(f1, wheretosave, 'epsc');

wheretosave = sprintf("graphs/f_lyap_%s", figname);
saveas(f2, wheretosave, 'epsc');