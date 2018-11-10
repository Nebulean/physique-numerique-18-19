clear
%% Définition des constantes
g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);

Omega = omega0;

n = 500; % On choisit le nombre de de pas de temps que l'on veut avoir dans chaque période.

tfin = 30000*2*pi/Omega; % Le cas où i = 10000 est celui de base, mais j'ai pris un peu plus.

dt = 2*pi/(n*Omega);

%% CHOIX DE DIFFERENTES CONDITIONS INITIALES
% Big brother is watching you. -> A GARDER: CAS NON-CHAOTIQUE (?)
% theta0 = 0.;
% thetadot0 = 1e-2;
% figname = "bigbrother";

% Une donut mince et compacte
% theta0 = pi/2;
% thetadot0 = 1e-2;

% Un donut moins mince avec des motifs plus joli.
% theta0 = pi/3;
% thetadot0 = pi;

% Motifs complétement chaotique. -> A GARDER: CAS CHAOTIQUE
% theta0 = pi;
% thetadot0 = 1e-2;
% figname = "chaos";

% Big fat donut très joli. -> A GARDER: CAS NON-CHAOTIQUE
% theta0 = 1e-6;
% thetadot0 = pi;
% figname = "fatdonut";

% Petits mouvements -> A GARDER: CAS NON-CHAOTIQUE
theta0 = 1e-6;
thetadot0 = 0.;
figname = "littlemoves";

%% On lance la simulation
cmd = sprintf("./Exercice3 configuration.in Omega=%s d=0.04 kappa=0. theta0=%s thetadot0=%s dt=%s tFin=%s sampling=%s output=e_poincare.out", num2str(Omega), num2str(theta0), num2str(thetadot0), num2str(dt), num2str(tfin), num2str(n))
system(cmd);



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