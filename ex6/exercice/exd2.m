clear all;
format long;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);

%% Initialisation
valb=0.02;
valR=0.12;

N=100;
valN1=N;
valN2=2*N;

%% Simulation
param = sprintf("output=outputdii trivial=false b=%f R=%f N1=%i N2=%i",valb,valR,valN1,valN2);
cmd = sprintf("./Exercice6 configuration.in %s",param);
disp(cmd);
system(cmd);

%% On importe et traite les données
data = load("outputdii_rholib_divEr_divDr.out");
rmidmid = data(:,1);
rholib = data(:,2);
divDr = data(:,4);

%% Graphs
figdiff=figure;
hold on;

grid on;
box on;

plot(rmidmid, divDr, 'x');
plot(rmidmid, rholib, '+');

zoomOfPlot(figdiff, 0.5, 0.3, 0.3, 0.3, [0.015, 0.025], [0, 2000])

xlabel("$r~[m]$");
ylabel("$\rho/\epsilon_0~[V/m^2]$");

l=legend(["$\nabla\cdot \mathbf{D}/\epsilon_0$", "$\rho_{lib}/\epsilon_0$"]);

set(gca, 'LineWidth',1.5);
set(gca, 'fontsize',25);

hold off;


saveas(figdiff, "graphs/exdii-diff", "epsc");


%% Quelques fonctions annexes
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
    main = get(gca, 'Position')
    g=copyobj(gca,fig)
    set(g, 'Position', [originx, originy, lengthx, lengthy]);
    set(g, 'XLim', limx);
    set(g, 'YLim', limy);
    set(g, 'LineWidth', 1.5);
end



function [polynome, slope] = poly_approx(x, y, ordre, steps)
    pf = polyfit(x, y, ordre);
    slope = pf(1);
    T = linspace(min(x), max(x), steps);

    n = ordre + 1;

    polynome = zeros(2,length(T));
    for i=1:n
       polynome(2,:) = polynome(2,:) + pf(i)*T.^(n-i);
    end

    polynome(1,:) = T;
    
    polynome=polynome.'; % transposition pour améliorer l'utilisation.
end