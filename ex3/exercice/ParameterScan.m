% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 20; % Nombre de simulations a faire

dt = ones(1,nsimul); % TODO: Choisir des valeurs de dt pour faire une etude de convergence
Omega = ones(1,nsimul); % TODO: Choisir des valeurs de Omega pour trouver la resonance

paramstr = 'dt'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = dt; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

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

if strcmp(paramstr, 'dt')
    error = zeros(1,nsimul);
elseif strcmp(paramstr, 'Omega')
    Emax = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    
    if strcmp(paramstr, 'dt')
        error(i) = 0; % TODO: Calculer l'erreur a partir de l'output
    elseif strcmp(paramstr, 'Omega')
        Emax(i) = 0; % TODO: Calculer le maximum de l'energie
    end
end

%% Figures %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    figure
    loglog(dt, error, 'k+')
    xlabel('\Delta t')
    ylabel('Erreur sur \theta(t_{fin}) [rad]')
    grid on
elseif strcmp(paramstr, 'Omega')
    figure
    plot(Omega, Emax, 'k-+')
    xlabel('\Omega [rad/s]')
    ylabel('max(E_{mec}(t)) [J]')
    grid on
end


