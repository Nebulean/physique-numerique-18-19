clear
%% Définition des constantes
g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);

Omega = 2*omega0;

n = 500; % On choisit le nombre de de pas de temps que l'on veut avoir dans chaque période.

tfin = 1000*2*pi/Omega;

dt = 2*pi/(n*Omega);

%% CHOIX DE DIFFERENTES CONDITIONS INITIALES
% MK8 % -> A GARDER PETIT MOUVEMENTS
theta01 = 0.;
thetadot01 = 1e-2;
% figname = "MK8";

% theta0 = 3*pi/4;
% thetadot0 = 0;

theta02 = pi/3;
thetadot02 = pi;

% Motifs complétement chaotique. -> A GARDER: CAS CHAOTIQUE
theta03 = pi;
thetadot03 = 1e-2;
% figname = "chaos1";

%-> A GARDER: CAS NON-CHAOTIQUE
% theta0 = 1e-6;
% thetadot0 = 1.;
% figname = "peanut";

% Petits mouvements -> A GARDER: CAS NON-CHAOTIQUE
theta04 = 1e-6;
thetadot04 = 0.;
% figname = "lemon";

% theta0 = 5*pi/6;
% thetadot0 = 0.;
% figname = "chaos2";

%% On lance la simulation
cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.05 kappa=0.1 theta0=%0.15f thetadot0=%0.15f dt=%0.15f tFin=%f sampling=%d output=f_poincare1.out", Omega, theta01, thetadot01, dt, tfin, n)
system(cmd);
disp(cmd);
cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.05 kappa=0.1 theta0=%0.15f thetadot0=%0.15f dt=%0.15f tFin=%f sampling=%d output=f_poincare2.out", Omega, theta02, thetadot02, dt, tfin, n)
system(cmd);
disp(cmd);
cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.05 kappa=0.1 theta0=%0.15f thetadot0=%0.15f dt=%0.15f tFin=%f sampling=%d output=f_poincare3.out", Omega, theta03, thetadot03, dt, tfin, n)
system(cmd);
disp(cmd);
cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.05 kappa=0.1 theta0=%0.15f thetadot0=%0.15f dt=%0.15f tFin=%f sampling=%d output=f_poincare4.out", Omega, theta04, thetadot04, dt, tfin, n)
system(cmd);
disp(cmd);

%% On load les différentes valeurs.
d1 = load("f_poincare1.out");

theta1 = d1(:,2);
thetadot1 = d1(:,3);

d2 = load("f_poincare2.out");

theta2 = d2(:,2);
thetadot2 = d2(:,3);

d3 = load("f_poincare3.out");

theta3 = d3(:,2);
thetadot3 = d3(:,3);

d4 = load("f_poincare4.out");

theta4 = d4(:,2);
thetadot4 = d4(:,3);

%% On dessine la section de poincaré.
fig=figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);
theta1=wrapToPi(theta1);
theta2=wrapToPi(theta2);
theta3=wrapToPi(theta3);
theta4=wrapToPi(theta4);
plot(theta1, thetadot1, '.');
plot(theta2, thetadot2, '.');
plot(theta3, thetadot3, '.');
plot(theta4, thetadot4, '.');

xlabel("$\theta$ [rad]");
ylabel("$\dot{\theta}$ [rad/s]");
lgd=legend("$\theta_0=0;\dot{\theta}_0=10^{-2}$", "$\theta_0=\pi/3;\dot{\theta}_0=\pi$", "$\theta_0=\pi;\dot{\theta}_0=10^{-2}$", "$\theta_0=10^{-6};\dot{\theta}_0=0$");
lgd.Location='northwest'
grid on;

hold off;

wheretosave = sprintf("graphs/f_poincare");
% saveas(fig, wheretosave, 'epsc');
print(fig, wheretosave,'-dpng','-r600');