%clear
%% Configuration des simulations
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);
Omega = omega0;

tfin = 20*(2*pi/Omega);

nsimul = 100; % Nombre de simulations a faire
nsteps = logspace(2,5,nsimul); % Nombre d'iterations entier de 10^2 a 10^4

dt = tfin ./ nsteps;
paramstr = 'dt';
param=dt;

%% SIMULATIONS
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(dt(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s Omega=%s d=0.04 kappa=0. theta0=0. thetadot0=1e-2 tfin=%s g=9.81', repertoire, executable, input, paramstr, dt(i), output{i}, num2str(Omega), num2str(tfin));
    disp(cmd) % Affiche la commande exécutée
    system(cmd); % Execute les simulations
end


%% CALCULS
thetafin = zeros(1,nsimul);
for i=1:nsimul
    output = sprintf('dt=%s.out', num2str(dt(i)));
    d = load(output);
    
    thetafin(i) = d(end,2);
end


%% Graph de la figure
f = figure
hold on;

p = plot(dt, thetafin, '.','LineWidth',1.2);

xlabel('$\Delta t$ [s]');
ylabel('Final angle $\theta$ [rad]')

set(gca, 'XScale','log');
%set(gca, 'YScale','log');
set(gca, 'fontsize', 20);
grid on;

%ymin=0.9e-3; ymax=1e-1;
%ylim([ymin ymax]);

hold off;

saveas(f, 'graphs/e_conv','epsc');