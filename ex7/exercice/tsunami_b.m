%% Préparation de(s) simulation(s)w
warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
schema = "B";
omega = 2*pi/900;
n_stride = 10;
Npoints = 1000;
tfin = 15000;
pulse = "true";
%Différents choix de xa.
nsimul = 6;
xa = round(linspace(250000, 350000, nsimul));
for i=1:nsimul
    param(i) = sprintf("schema=%s n_stride=%i Npoints=%i tfin=%i omega=%f pulse=%s xa=%i", schema, n_stride, Npoints, tfin, omega, pulse, xa(i));
    output(i) = sprintf("tsunami_b_%s_%i", schema, xa(i));
end


set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);


%% Simulations
for i=1:nsimul
    cmd = sprintf("./Exercice7 configuration_tsunami.in output=%s %s", output(i), param(i));
    disp(cmd);
    system(cmd);
end

x1wave = {};
t1wave = {};
u_calc = {};
ampl = {};
xampl = {};

for j=1:nsimul
    %% Chargement des données
    dataf = load(sprintf("%s_f.out", output(j)));
    datau = load(sprintf("%s_u.out", output(j)));
    dataE = load(sprintf("%s_E.out", output(j)));

    %% Traitement des données
    f = dataf(:,2:end);
    t = dataf(:,1);

    x = datau(:,1);
    
    % On veut suivre une unique vague, du coup, je dois trouver le maximum de
    % la première vague. On selectionne le dernier pic pour chaque instant.
    x1wave{j} = [];
    t1wave{j} = [];
    ampl{j} = [];
    xampl{j} = [];
    tE{j} = [];
    E{j} = [];
    
    
    isFinish = false;
    oldMax = -1.0;
    for i=1:length(f(:,1))
        % On trouve l'ensemble des pics de la fonction à un instant.
        [pks, indexes] = findpeaks(f(i,:));
        
        
        % On ne considère que les cas où il existe des pics (où c'est pas
        % nul), et les cas où la première vague n'a pas atteint la fin.
        if ~isempty(indexes) && ~isFinish
            % On garde que les trois points autour du premier pic.
            index = [indexes(end)-2, indexes(end)-1, indexes(end), indexes(end)+1, indexes(end)+2];
            
            % On interpol autour du résultat de ces trois points, si ils
            % existent.
            % On ne considère pas les bords.
            if ifIndexExists(x, index)
                ordre = 2;
                nbPoints = 50;
                % On trouve les informations du fit
                fit = polyfit(x(index), f(i,index).', ordre);
                
                % On créé un ensemble de points.
                T = linspace(min(x(index)), max(x(index)), nbPoints);
                
                % On "simule" la fonction que ça doit donner.
                interpol = zeros(2, length(T));
                for k=1:ordre+1
                   interpol(2,:) = interpol(2,:) + fit(k)*T.^(ordre+1-k);
                end
                interpol(1,:) = T;
                
                % Maintenant qu'on a l'interpolation d'ordre 2, on va
                % trouver le maximum.
                [M, idx] = max(interpol(2,:));
                ampl{j}(i) = M;
                xampl{j}(i) = interpol(1,idx);
                
                % On check la pente des X derniers points pour savoir si la
                % première vague est passée.
                XLastPoints = 3;
                [poly, slope] = poly_approx(x(end+1-XLastPoints:end), f(i,end+1-XLastPoints:end).', 1, 2, false);
                
                if slope > 0
                    isFinish = true;
                end
                % Si la vague atteint la fin, on stop pour cette
                % simulation.
                if ~isFinish
                    x1wave{j}(i) = interpol(1, idx); 
                    t1wave{j}(i) = t(i);
                    tE{j}(i) = dataE(i,1);
                    E{j}(i) = dataE(i,2);
                end
            end
        end
    end
    
    u_calc{j} = diff(x1wave{j})./diff(t1wave{j});
    x1wave{j} = x1wave{j}(1:end-1);
end


% On créé la valeur théorique
g=9.81;
uth = sqrt(g*h(x, xa));

%% On fait des graphes
fig_u=figure
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);


%plot(x, uth, '-', 'LineWidth', 1.5); % Vitesse théorique
for j=1:nsimul
    plot(x1wave{j}(1:end), u_calc{j}, ':', 'linewidth', 2); % Vitesse expérimentale
end

xlabel("$x~[m]$");
ylabel("$u~[m/s]$");

lgd = {};
for j=1:nsimul
   lgd{j} = sprintf("$x_a=%i km$", xa(j)/1000 ); 
end

legend(lgd, 'location', 'best');

grid on;
box on;

hold off;


%% Fig amp B
fig_f_B = figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

fA = (g*h(xampl{1}, xa)).^(-1/4);

% Choix de la constante en fonction de l'expérience:
%A0 = 1/0.05975; % Pour xa = 100
%A0 = 1/0.05975; % Pour xa = 250000
A0 = 1/0.05975; % Pour xa = 300000
%plot(xampl{1}, A0.*fA, '-', 'linewidth', 2)
%plot(xampl{1}, ampl{1} , '-', 'linewidth', 2);

for j=1:nsimul
    plot(xampl{j}, ampl{j} , ':', 'linewidth', 2);
end

lgd = {};
for j=1:nsimul
   lgd{j} = sprintf("$x_a=%i km$", xa(j)/1000 ); 
end

legend(lgd, 'location', 'best');

ylim([-0.5, 4.5]);

xlabel("$x~[m]$");
ylabel("$f~[m]$");
box on;
grid on;

hold off;

%% Quelques fonctions
function res = diff3pts(vec, order)
    %res = zeros(length(vec) - 2);

    for i=1:length(vec)-order
        res(i) = vec(i+order) - vec(i);
    end
end

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
    %   1, il est suffisant de prendre 2 points.
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


function res = ifIndexExists(vec, index)
    res = true;
    for i=1:length(index)
        %fprintf("current index is %i\n", index(i))
        if (length(vec) < index(i)) | (index(i) < 1)
            res = false;
            fprintf("Index %i does not exist.\n", index(i));
            break;
        end
    end
end

function res = h(x, xa)
    xb = 370000;
    xc = 430000;
    xd = 600000;
    L = 800000;
    
    hrecif = 20;
    hocean = 8000;
    
    res = [];
    
    for i=1:length(x)
        if (0 <= x(i)) & (x(i) < xa)
            res(i) = hocean;
        elseif (xa <= x(i)) & (x(i) < xb)
            res(i) = hocean + (hrecif-hocean) * sin(pi*(x(i)-xa)/(2*(xb-xa)))^2;
        elseif (xb <= x(i)) & (x(i) < xc)
            res(i) = hrecif;
        elseif (xc <= x(i)) & (x(i) < xd)
            res(i) = hrecif - (hrecif - hocean) * sin(pi*(xc-x(i))/(2*(xc-xd)))^2;
        elseif (xd <= x(i)) & (x(i) <= L)
            res(i) = hocean;
        end
    end
end



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