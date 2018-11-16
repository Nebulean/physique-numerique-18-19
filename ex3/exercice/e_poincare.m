clear
%% Définition des constantes
g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);

Omega = omega0;

n = 100; % On choisit le nombre de de pas de temps que l'on veut avoir dans chaque période.

tfin = 10000*2*pi/Omega; % Le cas où i = 10000 est celui de base, mais j'ai pris un peu plus.

dt = 2*pi/(n*Omega);

%% CHOIX DE DIFFERENTES CONDITIONS INITIALES
% Big brother is watching you. -> A GARDER: CAS NON-CHAOTIQUE (?)
theta0 = 0.;
thetadot0 = 1e-2;
figname = "bigbrother";

% Une donut mince et compacte
% theta0 = pi/2;
% thetadot0 = 1e-2;

% Un donut moins mince avec des motifs plus joli.
% theta0 = pi/3;
% thetadot0 = pi;

% Motifs complétement chaotique. -> A GARDER: CAS CHAOTIQUE
% theta0 = 3.12;
% thetadot0 = 1e-2;
% figname = "chaos";

% Big fat donut très joli. -> A GARDER: CAS NON-CHAOTIQUE
% theta0 = 1e-6;
% thetadot0 = pi;
% figname = "fatdonut";

% Petits mouvements -> A GARDER: CAS NON-CHAOTIQUE
% theta0 = 1e-6;
% thetadot0 = 0.;
% figname = "littlemoves";

%% On lance la simulation
cmd = sprintf("./Exercice3 configuration.in Omega=%0.15f d=0.04 kappa=0. theta0=%f thetadot0=%f dt=%0.15f tFin=%f sampling=%d output=e_poincare.out", Omega, theta0, thetadot0, dt, tfin, n)
system(cmd);
disp(cmd);



%% On load les différentes valeurs.
d = load("e_poincare.out");

theta = d(:,2);
thetadot = d(:,3);

%% On dessine la section de poincaré.
fig=figure;
hold on;

plot(thetadot, theta, '.');

xlabel("$\dot{\theta}$ [rad/s]");
ylabel("$\theta$ [rad]")
grid on;

hold off;

wheretosave = sprintf("graphs/e_poincare_%s", figname);
saveas(fig, wheretosave, 'epsc');