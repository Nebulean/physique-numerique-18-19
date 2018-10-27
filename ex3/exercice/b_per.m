%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 1000; % Nombre de simulations a faire

theta0 = linspace(0.0001,3.1415,nsimul);

paramstr = 'theta0';
param = theta0;

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

%Copie des paramètres de configuration.in
tfin = 20;
dt = 0.02;
g = 9.81;
L = 0.1;

omega0 = sqrt(g/L);

period = zeros(nsimul,1); % vecteur qui rassemblera les périodes de chaque simulation
periodth = zeros(nsimul,1); % valeurs théoriques

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    t = data(:,1);
    theta = data(:,2);
    
    times = []; % création du vecteur des t auxquels theta est maximum
    
    for j = 2:length(t)-1
       if (theta(j)>theta(j-1) && theta(j)>theta(j+1))
           times = [times; t(j)]; % on ajoute le temps à chaque maximum de theta
       end
    end
    
    times = flipud(times); % on inverse le sens du vecteur pour itérer dessus
    
    for j = 1:length(times)-1
        times(j) = times(j) - times (j+1); % on calcule chaque période en faisant la différence entre les temps à chaque 2 maxima consécutifs
    end
    
    %disp(length(times));
    %disp(times);
    
    period(i,1) = sum(times)/length(times); % on rassemble la moyenne des périodes de chaque simulation
    periodth(i,1) = 4/omega0 * ellipticK(sin(theta0(i)/2)*sin(theta0(i)/2)); % valeur théorique de chaque simulation
end

%disp(period);

%% Figures %%
%%%%%%%%%%%%%

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

figure
p1 = plot(theta0, period, theta0, periodth);
set(p1, 'LineWidth', 1.5);
set(gca,'fontsize',20);
xlabel('$\theta_0$ [rad]')
ylabel('Period T [s]')
grid on
legend('Verlet','Theoretical')