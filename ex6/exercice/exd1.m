clear all;
format long;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);

%% Initialisation
valb=0.02;
valR=0.12;

% nsimul = 10; % Pour convergence
% N=logspace(1,3,nsimul); % Pour convergence
nsimul = 4; % Pour evolution
N = linspace(10,100,nsimul) % Pour evolution

N=int32(N);

output = {};
for i = 1:nsimul
    output{i} = sprintf("output%i",i);
end

%% On prépare et lance les simulations.
for i=1:nsimul
    valN1=N(i);
    valN2=2*N(i);
    
    param = sprintf("output=%s trivial=false b=%f R=%f N1=%i N2=%i",output{i},valb,valR,valN1,valN2);
    cmd = sprintf("./Exercice6 configuration.in %s",param);
    disp(cmd);
    system(cmd);
end


%% On importe et traite les données
for n=1:nsimul
    clear data index;
    data = load(sprintf("%s_phi.out", output{n}));
    
    % On cherche l'index où se trouve b.
    index = 0;
    for i=1:length(data)
       if(data(i,1) == valb)
          index = i;
       end
    end
    
    phi(n) = data(index,2);
    
    data = load(sprintf("%s_Er_Dr.out", output{n}));
    Er(n) = data(index,2); % On prend que la valeur en b
    rmid{n} = data(:,1); % On prend tout
    Ertot{n} = data(:,2); % On prend tout
    
    data = load(sprintf('%s_phi.out', output{n}));
    r{n} = data(:,1); % On prend tout
    phitot{n} = data(:,2); % On prend tout
end


%% Graphs
figphi=figure;
hold on;

set(gca, 'fontsize', 22);
set(gca, 'LineWidth',1.5);

plot((1./single(N)).^2,phi,'+', 'MarkerSize',10);
[pn, slope] = poly_approx((1./single(N)).^2,phi,1,2);

plot(pn(:,1), pn(:,2), 'LineWidth',1.2);

xlabel("$1/N^2~[-]$");
ylabel("$\phi(r)~[V]$");

legend(["Data", "Linear fit"])

grid on;
box on;

hold off;


% Pas faire plus que 4.
figEvoEr=figure;
hold on;

grid on;
box on;

for i=1:nsimul
    plot(rmid{i}, Ertot{i}, '.-.');
    lgd{i} = sprintf("N1=%i, N2=%i", N(i), 2*N(i));
end

zoomOfPlot(figEvoEr, 0.4, 0.27, 0.3, 0.3, [0.012, 0.021], [48, 55])

set(gca, 'fontsize', 22);
set(gca, 'LineWidth',1.5);

legend(lgd);

xlabel("$r~[m]$");
ylabel("$E(r)~[V/m]$")

hold off;




figEvoPhi=figure;
hold on;

grid on;
box on;

for i=1:nsimul
    plot(r{i}, phitot{i}, '.-.');%, 'LineWidth',1);
    lgd{i} = sprintf("N1=%i, N2=%i", N(i), 2*N(i));
end

zoomOfPlot(figEvoPhi, 0.55, 0.35, 0.25, 0.25, [0.015, 0.025], [0.3,0.4]);

set(gca, 'fontsize', 22);
set(gca, 'LineWidth',1.5);

legend(lgd);

xlabel("$r~[m]$");
ylabel("$\phi(r)~[V]$")


hold off;



saveas(figphi, "graphs/exd1-convPhi","epsc");
saveas(figEvoPhi, "graphs/exd1-evoPhi","epsc");
saveas(figEvoEr, "graphs/exd1-evoEr", "epsc");


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
    set(g, 'LineWidth', 1);
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