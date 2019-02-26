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

dt = logspace(-5,-3, nsimul);

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
    yp = .02;
    Tp = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'dt'))
        data = load([output{i} '_T.out']);
        N = sqrt(length(data));
        Y = data(:,2);
        X = data(:,1);
        T = data(:,3);
        
%         Xap =  min(abs(X-xp))
%         if (Xap <= xp)
%             Xlow = Xap
%         else
%             Xlow = Xap-1
%         end
%         Xhigh = Xlow+1;
%         
%         Yap =  min(abs(Y-xp));
%         if (Yap <= yp)
%             Ylow = Yap
%         else
%             Ylow = Xap-1
%         end
%         Yhigh = Ylow+1;

%         Xlow = floor(xp*N/L)
%         Xhigh = Xlow+1
%         Ylow = floor(yp*N/L)
%         Yhigh = Ylow+1

        Xid = xp*(N-1)/L+1
        Yid = yp*(N-1)/L+1

%         
%         T1 = T(Xlow+N*Ylow+1)
%         T2 = T(Xlow+N*Yhigh+1)
%         T3 = T(Xhigh+N*Ylow+1)
%         T4 = T(Xhigh+N*Yhigh+1)
%         
        Tid = (Xid-1)*N+Yid
        Tp(i) = T(Tid)
    end
end

%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paramstr,'dt'))
    figure
    plot(dt,Tp,'k+')
    xlabel('\Delta t [s]')
    ylabel(sprintf('T(%0.2f,%0.2f) [°C]',xp,yp))
    grid on
end



