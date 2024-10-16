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
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 20; % Nombre de simulations a faire

nsteps = round(logspace(2,4,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4
tfin = 1.09321e-7;
dt = tfin ./ nsteps;

paramstr = 'nsteps'; % Nom du parametre a scanner
param = nsteps; % Valeurs du parametre a scanner

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
% Copie direct des valeurs de configuration.in
tfin   = 1.09321e-7;
q      = 1.6022e-19;
m      = 1.6726e-27;
B0     = 3;
E      = 0;
Kappa  = 0;
vx0    = 0;
vy0    = 4e5;
x0     = -vy0*m/(q*B0);
y0     = 0;

v0 = vy0;
omega=q*B0/m;

error = zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    t = data(:,1);
    x = data(:,2);
    y = data(:,3);
    x_th = v0*sin(omega*t)/omega; % TODO: Entrer la vraie solution analytique en fonction du temps
    y_th = -v0*cos(omega*t)/omega; % TODO: Entrer la vraie solution analytique en fonction du temps
    error(i) = max(sqrt((x-x_th).^2+(y-y_th).^2));
end

figure
loglog(dt, error, 'k+')
xlabel('\Delta t')
ylabel('Maximum de l''erreur sur la position')
grid on
