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
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

% low=10;
% high=100;
nsimul=60;

Npoints = round(logspace(1,3, nsimul));

paramstr = 'Npoints'; % Nom du parametre a scanner
param = Npoints; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g cb_droit=sortie tfin=2. output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd);
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%
omega=5.;
tp=1.5;
xp=5.;
L=20.;
u=6.;


f=zeros(1,nsimul);
f_th=-sin((omega/u)*xp-omega*tp);
time=zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'Npoints'))
        dataf = load([output{i} '_f.out']);
        datau = load([output{i} '_u.out']);
        time = dataf(:,1);
        pos = datau(:,1);
%         [T,X]=meshgrid(time,pos);
        func = dataf(:,2:end)
%         [diff,tid] = min(abs(time-tp));
        
%         dx=L/(Npoints(i)-1);
%         xid=round(xp/dx)
        
        f(i)=griddata(pos,time,func,xp,tp,'linear')
%         f(i)=interp2(pos,time,func,xp,tp,'spline')
    end
end

err = abs(f-f_th);

%% Figures %%
%%%%%%%%%%%%%


if(strcmp(paramstr,'Npoints'))
    f1=figure;
    
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 18);
    set(gca, 'fontsize', 25);
    set(gca, 'LineWidth',1.5);
    
    hold on
    plot(Npoints,err,'k+');
    
    %resulto pimpagu no jutsu
    errfit=err;
    Nfit=Npoints
    for i=1:length(err);
        if err(i)<1e-10
            errfit(i)=[];
            Nfit(i)=[];
        end
    end
    
    [fit, slope] = poly_approx(Nfit, errfit, 1, 2, true);
    plot(fit(:,1), fit(:,2), '-');
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlabel('$N$')
    ylabel('$|f-f_{th}|$ [m]')
	legend(["Data", "slope ="+num2str(slope)],'Location','southeast');
    box on;
    grid on
    hold off;
end

saveas(f1, "graphs/convN","epsc");

%% Fonction

function [polynome, slope] = poly_approx(x, y, ordre, steps, isLog)
    % =================== POLY_APPROX ===================================
    % RESUMÉ: Permet de faire une approximation d'ordre n d'un set de
    % données.
    %
    % USAGE: Il suffit d'appeller la fonction.
    %
    % PARAMETRES:
    %   - (x/y): Set de donnée de même longueur.
    %   - order: L'ordre du fit.
    %          Valeurs acceptables: Tous les entiers plus grand que 0.
    %   - steps: Nombre de points à sortir. Pour une approximation d'ordre
    %   1, il est suffisent de prendre 2 points.
    %          Valeurs acceptables: Tout nombre entier supérieur à 2.
    %   - isLog: Switch permettant de choisir entre un graphe linéaire et un
    %   graphe log. 'true' si graphe log, 'false' si graphe linéaire.
    % ==================================================================
    if isLog == true
        x=log10(x);
        y=log10(y);
    end
    
    % Compute the fit
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
    
    if isLog == true
       polynome = 10.^polynome; 
    end
end