%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
delta = 0;
x0 = 0;
sigma_norm = 0.06;
n = 14;
tfin = 5000;
Ninters = 300;

dt = 1;

cmd = sprintf("./Exercice8 configuration.in output=i_ptot xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", xL, xR, omega, delta, x0, sigma_norm, n, tfin, Ninters, dt);


%% Simulations
disp(cmd);
system(cmd);


%% Traitement des données.

psi2 = load("i_ptot_psi2.out");

data = load("i_ptot_obs.out")
t = data(:,1);
probG = data(:,2);
probD = data(:,3);
E = data(:,4);
xmoy = data(:,5);
x2moy = data(:,6);
pmoy = data(:,7);
p2moy = data(:,8);
delx = data(:,9);
delp = data(:,10);

data = load("i_ptot_pot.out");
x = data(:,1);
V = data(:,2);

[X, T] = meshgrid(x,t);


%% Figures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);


figPtot = figure;
hold on;
%dt = 5, bonne droite.
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, probD + probG, '-', 'linewidth', 2);

xlabel("$t~[t_P]$");
ylabel("$\sum P(t)$");

% Old version
% plot(t, probD + probG, '-', 'linewidth', 5);
% test = [probG, probD];
% area(t, test);
% 
% xlabel("$t~[t_P]$");
% ylabel("$P(t)$");
% 
% legend("$\Sigma P(t)$", "$P_{x<0}$", "$P_{x>0}$")
% 
% ylim([0, 1.5]);

grid on;
box on;

hold off;

%%%%%%%%%%%%%%%%%%%%%%%

figE = figure;
hold on;
%dt = 1, belle courbe.
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, E./E(1), '-', 'linewidth', 2);

xlabel("$t~[t_P]$");
ylabel("$E(t)/E(0)~[-]$");

grid on;
box on;

hold off;

%%%%%%%%%%%%%%%%%%%%%%%

figHeisenberg = figure;
hold on;
%dt = X, bonne droite.
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, delx.*delp, '-', 'linewidth', 2);
line([0, 5000], [1/2 1/2], 'linewidth', 2, 'color', [0.8, 0.2, 0.2]);
grid on;
box on;

% zoomOfPlot(figHeisenberg, 0.3, 0.7, 0.18, 0.18, [501, 549], [0.4991, 0.5009]);

xlabel("$t~[t_P]$");
ylabel("$\langle \Delta x \rangle(t) \langle \Delta p \rangle(t)~[\ell_P~m_P~c]$");

legend(["$\langle \Delta x \rangle(t) \langle \Delta p \rangle(t)$", "$\hbar/2$"])

ylim([0.48, 0.95]);

hold off;


%%%%%%%%%%%%

figEvo=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

pcolor(X, T, psi2);
colbar=colorbar;
shading interp;

box on

xlabel("$x~[\ell_P]$");
ylabel("$t~[t_P]$");
ylabel(colbar, "$|\psi(x,t)|^2$", 'interpreter', 'latex', 'fontsize', 25);

hold off;



%% Fonctions
function zoomOfPlot(fig, originx, originy, lengthx, lengthy, limx, limy)
    % ================== ZOOMOFPLOT ====================================
    % RESUMÉ: Permet de faire un zoom sur une zone du graphe actuel.
    %
    % USAGE: Il suffit d'appeller la fonction lorsqu'on est dans une figure
    % (avant le hold off du coup)
    %
    % PARAMETRES:
    %   - fig: nom de la figure courante. Cette figure doit donc avoir un
    %          nom (f=figure par exemple).
    %   - origin(x/y): origine (bord bas-gauche) de la figure zoomée.
    %          Valeurs acceptables: [0 à 1].
    %   - length(x/y): longueur de la figure zoomée.
    %          Valeurs acceptables: [0 à 1].
    %   - lim(x/y): zone d'intéret pour le zoom.
    %          Valeurs acceptables: limites des axes.
    %
    % REMARQUE: Cette fonction fait une copie de tout ce qui ce qui a été
    % fait dans la figure jusqu'à l'appel. Il faut donc mettre les label,
    % legendes et autres après l'appel à cette fonction.
    % ==================================================================
    % On fait une boite
    lx = abs(limx(2) - limx(1));
    ly = abs(limy(2) - limy(1));
    rectangle('Position', [limx(1), limy(1), lx, ly])
    
    % On fait le plot zoomé.
    g=copyobj(gca,fig)
    set(g, 'Position', [originx, originy, lengthx, lengthy]);
    set(g, 'XLim', limx);
    set(g, 'YLim', limy);
    set(g, 'LineWidth', 1);
end