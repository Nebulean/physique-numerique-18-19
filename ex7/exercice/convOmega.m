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

low=0.;
high=6;
nsimul=600;

omega = linspace(low,high, nsimul);

paramstr = 'omega'; % Nom du parametre a scanner
param = omega; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g cb_gauche=harmonique cb_droit=fixe tfin=400 n_stride=10 ecrire_f=false output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd);
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%
u = 6.
L = 20.

maxE=zeros(1,nsimul);
time=zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'omega'))
        data = load([output{i} '_E.out']);
        time = data(:,1);
        E = data(:,2);
        
        maxE(i)=max(E);
        
    end
end

%% Figures %%
%%%%%%%%%%%%%

u = 6.
L = 20.

if(strcmp(paramstr,'omega'))
    f1=figure;
    
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 18);
    set(gca, 'fontsize', 25);
    set(gca, 'LineWidth',1.5);
    
    hold on
    plot(omega,maxE);

    for i=1:6
        omega0=u/L*i*pi;
        line([omega0 omega0],[0 15e4], 'color', 'red');
    end
    
    xlabel('Frequency of excitation [rad/s]')
    ylabel('$E_{max}$ [$m^3$]')
	legend("Simulations", "$\omega_n$");
    box on;
    grid on
    hold off;
    
%     f2=figure;
%     
%     set(gca, 'fontsize', 25);
%     set(gca, 'LineWidth',1.5);
%     hold on
%     plot(omega(1:end/9),maxE(1:end/9));
% 
% %     for i=1:3
%         omega0=u/L*pi;
%         line([omega0 omega0],[0 15e4], 'color', 'red');
% %     end
%     
%     xlabel('$\omega$')
%     ylabel('$E_{max}$ [$m^3$]')
%     legend("Simulations", "Analytical eigenfrequencies");
%     box on;
%     grid on
%     hold off;
end

saveas(f1, "graphs/Eom","epsc");
% saveas(f1, "graphs/Eomzoom","epsc");

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