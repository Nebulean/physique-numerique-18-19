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
executable = 'Exercice6'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

low=10;
high=50;
nsimul=high-low+1;

N1 = linspace(low,high, nsimul);

paramstr = 'N1'; % Nom du parametre a scanner
param = N1; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g N2=%.15g output=%s', repertoire, executable, input, paramstr, param(i), param(i), output{i});
    disp(cmd);
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

R=.12;
phi_0th=R*R/4;

phi_0=zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'N1'))
        data = load([output{i} '_phi.out']);
        phi_0(i) = data(1,2);
    end
end

err = abs(phi_0-phi_0th);

%% Figures %%
%%%%%%%%%%%%%


if(strcmp(paramstr,'N1'))
    f=figure;
    
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 18);
    set(gca, 'fontsize', 25);
    set(gca, 'LineWidth',1.5);
    
    hold on
    plot(N1,err,'k+');
%     [fit, slope] = poly_approx(N1, err, 1, 2, true);
%     plot(fit(:,1), fit(:,2), '-');
    %fit = 10.^fit;
%     set(gca, 'YScale', 'log');
%     set(gca, 'XScale', 'log');
    %poly_approx(dt, Tp, 1, 2, true);
%     [P,slope]=poly_approx(dt, Tp, 1, 2);
%     plot(P(1,:),P(2,:));
    xlabel('$N1$')
    ylabel('$|\phi(0)-\phi(0)_th|$ [V]')
	legend(["Data", "slope ="+num2str(slope)]);
    box on;
    grid on
    hold off;
end

saveas(f, "graphs/c_convphi","epsc");

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