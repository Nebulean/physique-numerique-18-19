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
executable = 'Exercice5'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 10; % Nombre de simulations a faire

% POUR N=40:
dt = logspace(-5,-3, nsimul);

% POUR N=80:
% dt = logspace(-5,-4, nsimul);

paramstr = 'dt'; % Nom du parametre a scanner
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

L = .1;

% Parcours des resultats de toutes les simulations

if(strcmp(paramstr,'dt'))
    xp = .05;
    yp = .05;
    Tp = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'dt'))
        data = load([output{i} '_T.out']);
        N = sqrt(length(data));
        Y = data(:,2);
        X = data(:,1);
        T = data(:,3);
        

%         Xid = xp*(N-1)/L+1
%         Yid = yp*(N-1)/L+1
% %         
%         Tid = (Xid-1)*N+Yid

        G = griddata(X,Y,T,xp,yp);
        Tp(i) = G
    end
end

%% Figures %%
%%%%%%%%%%%%%


if(strcmp(paramstr,'dt'))
    f=figure
    
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 18);
    set(gca, 'fontsize', 25);
    set(gca, 'LineWidth',1.5);
    
    hold on
    plot(dt,Tp,'k+');
    [P,slope]=poly_approx(dt, Tp, 1, 2);
    plot(P(1,:),P(2,:));
    xlabel('$\Delta t$ [s]')
    ylabel(sprintf('$T(%0.2f,%0.2f)$ [K]',xp,yp))
    legend("slope =" +num2str(slope));
    grid on
    hold off;
end

saveas(f, "graphs/b_conv40","epsc");
% saveas(f, "graphs/b_conv80","epsc");

%% Fonction

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
end
