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

%nsimul = 5; % Nombre de simulations: pour la figure 3.
nsimul = 50; % Nombre de simulations: pour les figures 1 et 2.

%nsteps = round(logspace(2.3,4,nsimul)); % pour la figure 3.
nsteps = round(logspace(2,6,nsimul)); % pour les figures 1 et 2.
tfin = 1.09321e-7;
dt = tfin ./ nsteps;

paramstr = 'nsteps'; % Nom du parametre a scanner
param = nsteps; % Valeurs du parametre a scanner

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
E      = 6e4;
Kappa  = 0;
vx0    = 0;
vy0    = 4e5;
x0     = -vy0*m/(q*B0);
y0     = 0;

v0 = vy0;
omega=q*B0/m;

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s schema=RK2 E=%s', repertoire, executable, input, paramstr, param(i), output{i}, num2str(E));
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

lastX = zeros(1,nsimul);
lastY = zeros(1,nsimul);
X={};
Y={};
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    t = data(:,3);
    x = data(:,2);
    y = data(:,3);
    lastX(i) = x(length(x));
    lastY(i) = y(length(y));
        
    X{i} = data(:,2);
    Y{i} = data(:,3);
end

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

f1 = figure;
hold on
p1 = plot(dt, lastX);
set(p1, 'LineWidth', 1.5);
xlabel('$\Delta$ t');
ylabel('$x(t_{final})$ [m]')
grid on
set(gca,'fontsize',20);
set(gca,'XScale','log');
%set(gca,'YScale','log');
hold off


f2 = figure;
hold on
p2 = plot(dt, lastY);
set(p2, 'LineWidth', 1.5);
xlabel('$\Delta$ t');
ylabel('$y(t_{final})$ [m]')
grid on
set(gca,'fontsize',20);
set(gca,'XScale','log');
%set(gca,'YScale','log');
hold off

f3 = figure
hold on
for i=1:nsimul
    p=plot(X{i}, Y{i},'LineWidth',1.4);
end
hold off

saveas(f1, 'graphs/app1_i_conv_x','epsc');
saveas(f2, 'graphs/app1_i_conv_y','epsc');
%saveas(f3, 'graphs/app1_i_conv_xy','epsc');