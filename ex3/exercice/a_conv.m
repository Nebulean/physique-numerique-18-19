%clear

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 20; % Nombre de simulations a faire
nsteps = logspace(3,5,nsimul); % Nombre d'iterations entier de 10^2 a 10^4
tfin = 20;
dt = tfin ./ nsteps;
paramstr = 'dt';
param=dt;

% SIMULATIONS
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(dt(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s kappa=0. d=0. theta0=1e-6 thetadot0=0. tfin=20 g=9.81', repertoire, executable, input, paramstr, dt(i), output{i});
    disp(cmd)
    system(cmd); % À commenter quand les simulations sont faites.
end


% CALCULS
% Solution analytique
A = 1e-6;
g = 9.81;
L = 0.1;
omega0 = sqrt(g/L);
%theta_th = A*cos(omega0*20);

% On calcul le tout
error = ones(1,nsimul);
thetafin = zeros(1,nsimul);
for i=1:nsimul
    output = sprintf('dt=%s.out', num2str(dt(i)));
    d = load(output);
    
    timefin = d(end,1);
    thetafin = d(end,2);
    thetath = A*cos(omega0*timefin);
    
    error(i) = abs(thetafin - thetath);
end

%% Graph de la figure
f = figure
hold on;

p = plot(dt, error, 'x','LineWidth',1.2);

xlabel('$\Delta t$ [s]');
ylabel('Error on the final angle $\theta$')

set(gca, 'XScale','log');
set(gca, 'YScale','log');
set(gca, 'fontsize', 20);
grid on;

hold off;

saveas(f, 'graphs/a_conv','epsc');