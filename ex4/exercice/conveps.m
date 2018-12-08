%% TASK
% Impémenter le pas de temps adaptatif dans le code et calculer la trajectoire. 
% Faire une étude de convergence, en variant ε, qui est la précision requise par pas de temps.
% Illuster comment le pas de temps ∆t change au cours de la simulation.
% Comparer le nombre de pas de temps nécessaires pour obtenir un résultat sur hmin de précision donnée
% avec le schéma à ∆t fixe.
%% BEFORE-SIMULATIONS
tFin = 172800; % 2 days in seconds.
dtad="true";
dt=60; % dt initial

nsimul = 20; % Number of simulations that we want. 8, because 172800 = 2^8 * 3^3 * 5^2.
% epsilon = zeros(nsimul, 1);
% for i=1:nsimul
%     epsilon(i,1) = 10/10^i;
%     %disp(sprintf("i = %d", i));
% end
epsilon = logspace(0, -5, nsimul)

%% SIMULATIONS
output = cell(1, nsimul)
for i=1:nsimul
     output{1, i} = sprintf("epsilon=%0.7f.out", epsilon(i));
     cmd = sprintf("./Exercice4 configuration.in epsilon=%0.15f dt=%0.15f dtad=%s tFin=%d output=%s", epsilon(i), dt, dtad, tFin, output{1,i});
     
     disp(cmd);
     system(cmd);
end

%% LOADING DATA
mindist = zeros(nsimul, 1);
nsteps = zeros(nsimul, 1);
RT = 6378.1 * 1000; % earth's radius
distTH = 10e3 + RT;
tfinal = zeros(nsimul, 1);
for i=1:nsimul
    data = load(output{1,i});

    Ex = data(:,3);
    Ey = data(:,4);
    Ax = data(:,11);
    Ay = data(:,12);
    
    t = data(:, 1);
    nsteps(i,1) = length(t);
    tfinal(i,1) = t(end);
    %% On calcul la distance minimale
    % Distance numérique avec les dt
    dist = sqrt((Ex - Ax).^2 + (Ey - Ay).^2);
    % Premièrement, on calcul le point où il y a la distance minimale.
    [tmp, index] = min(dist);
    if index+1 > length(dist)
        fit = polyfit(t(index-2:index), dist(index-2:index), 2);
    else
        fit = polyfit(t(index-1:index+1), dist(index-1:index+1), 2);
    end
    A = fit(1); B = fit(2); C = fit(3);
    
    % On calcul le minimum à l'aide de l'analyse
    mindist(i, 1) = abs(C - B^2/(4*A) - distTH); %TODO: Retirer les abs lorsque le résultat sera correcte
   t(1,i) = data(end, 1);
end

%% PLOTTING DATA
fig1=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(nsteps, mindist, 'x');
grid on;

% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');

xlabel("Number of steps [-]");
ylabel("$h_{min}$ [m]");

hold off;