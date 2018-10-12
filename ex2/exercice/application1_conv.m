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

schema = ["Euler", "EulerCromer", "RungeKutta2"];

% Quelques couleurs
darkblue = [0.1 0.1 0.8];
darkred = [0.8 0.1 0.1];
darkgreen = [0.2 0.8 0.2];
% le tout dans un array
color = {};
color{1} = darkblue;
color{2} = darkred;
color{3} = darkgreen;


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

f1 = figure;
hold on
for n=1:3
    %% Simulations %%
    %%%%%%%%%%%%%%%%%

    output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
    for i = 1:nsimul
        output{i} = [paramstr, '=', num2str(param(i)), '.out'];
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %s=%.15g output=%s schema=%s', repertoire, executable, input, paramstr, param(i), output{i}, schema{n});
        disp(cmd)
        system(cmd);
    end

    %% Analyse %%
    %%%%%%%%%%%%%
 

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

    
    p = plot(dt, error)
    set(p, 'LineWidth', 1.5);
    set(p, 'Color', color{n});
    xlabel('\Delta t')
    ylabel('Maximum de l''erreur sur la position')
    grid on
    set(gca,'fontsize',16);
end

legstr = ["Euler", "Euler Cromer", "Runge Kutta 2"];

l = legend(legstr);
set(l, 'Location', 'northwest');
set(l, 'FontSize', 14);
hold off



%% Vu que le premier ne donne pas beaucoup d'informations, on dessine un deuxième graph, sans Euler.
f2 = figure;
hold on
for n=2:3
    output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
    for i = 1:nsimul
        output{i} = [paramstr, '=', num2str(param(i)), '.out'];
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %s=%.15g output=%s schema=%s', repertoire, executable, input, paramstr, param(i), output{i}, schema{n});
        disp(cmd)
        system(cmd);
    end
    
    
    %% Analyse %%
    %%%%%%%%%%%%%
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

    
    p = plot(dt, error)
    set(p, 'LineWidth', 1.5);
    set(p, 'Color', color{n});
    xlabel('\Delta t')
    ylabel('Maximum de l''erreur sur la position')
    grid on
    set(gca,'fontsize',16);
end

legstr2 = ["Euler Cromer", "Runge Kutta 2"];

l2 = legend(legstr2);
set(l2, 'Location', 'northwest');
set(l2, 'FontSize', 14);
hold off

saveas(f1, 'graphs/app1_conv_ALL', 'epsc');
saveas(f2, 'graphs/app1_conv_noEuler', 'epsc');