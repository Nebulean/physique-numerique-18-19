

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice5'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

L = .1;
xa = .02;
xb = .04;
xc = .075;
xd = .085;
N = 40;
h = L/N

nsimul = int16((xc-xb)/(2*h)) % Nombre de simulations a faire

dx = linspace(0,(xc-xb)/2, nsimul);

paramstr = 'dx';
param = dx; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i),6)];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s xa=%.15g xb=%.15g xc=%.15g xd=%.15g output=%s', repertoire, executable, input, xa+dx(i), xb+dx(i), xc-dx(i), xd-dx(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

if(strcmp(paramstr,'dx'))
    xp = xb+(xc-xb)/2;
    yp = .05;
    Tp = zeros(1,nsimul);
    Fp = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'dx'))
        data = load([output{i} '_T.out']);
        N = sqrt(length(data));
        Y = data(:,2);
        X = data(:,1);
        T = data(:,3)-273.15;

        G = griddata(X,Y,T,xp,yp);
        Tp(i) = G;
        
        dataF = load([output{i} '_F.out']);
        Fx = dataF(:,3);
        Fy = dataF(:,4);
        Xmid = dataF(:,1);
        Ymid = dataF(:,2);
        Fnorm = sqrt(Fx.^2+Fy.^2);
        
        Fp(i) = griddata(Xmid,Ymid,Fnorm,xp,yp);
    end
end

%% Figures %%
%%%%%%%%%%%%%


if(strcmp(paramstr,'dx'))
    f1=figure;
    
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 18);
    set(gca, 'fontsize', 25);
    set(gca, 'LineWidth',1.5);
    
    hold on
    plot(dx,Tp,'k+');
    % [P,slope]=poly_approx(dt, Tp, 1, 2);
%     plot(P(1,:),P(2,:));
    xlabel('$\Delta x$ [m]')
    ylabel(sprintf('$T(%0.2f,%0.2f)$ [C]',xp,yp))
    % legend("slope =" +num2str(slope));
    grid on
    hold off;
    
    f2=figure;
    
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 18);
    set(gca, 'fontsize', 25);
    set(gca, 'LineWidth',1.5);
    
    hold on
    plot(1./sqrt(xc-xb-2*dx),Fp,'k+');
    [P,slope]=poly_approx(1./sqrt(xc-xb-2*dx), Fp, 1, 2)
    plot(P(1,:),P(2,:));
    xlabel('$1/\sqrt{d}$ [m$^{-1/2}$]')
    ylabel(sprintf('$|j|(%0.2f,%0.2f)$ [W/m$^2$]',xp,yp))
    legend("slope =" +num2str(slope));
    grid on
    hold off;
end

saveas(f1, "graphs/e_distT","epsc");
saveas(f2, "graphs/e_distF2","epsc");

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
