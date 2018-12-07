%% TASK
% Résoudre numériquement la trajectoire d’Apollo 13 avec un pas de temps ∆t fixe. [N.B. :
% Il suffit de simuler 2 jours]. Faire une étude de convergence de l’altitude minimale hmin
% et de vmax avec ∆t et comparer avec le résultat analytique.

%% BEFORE-SIMULATIONS
tFin = 172800; % 2 days in seconds.
dtad="true";
dt=60;

nsimul = 8; % Number of simulations that we want. 8, because 172800 = 2^8 * 3^3 * 5^2.
epsilon = zeros(nsimul, 1);
for i=1:nsimul
    epsilon(i,1) = 10/10^i;
    disp(sprintf("i = %d", i));
end
% epsilon = logspace(-7,-1,nsimul);

%% SIMULATIONS
output = cell(1, nsimul)
for i=1:nsimul
     output{1, i} = sprintf("epsilon=%0.7f.out", epsilon(i));
     cmd = sprintf("./Exercice4 configuration.in epsilon=%0.15f dt=%0.15f dtad=%s tFin=%d output=%s", epsilon(i), dt, dtad, tFin, output{1,i});
     
     disp(cmd);
     system(cmd);
end

%% LOADING DATA
fordist = zeros(nsimul, 7);
tfin = zeros(1, nsimul);
% forvmax = zeros(nsimul, 6);
for i=1:nsimul
   data = load(output{1,i});
   
   fordist(i, 1) = data(end,2); % dt
   fordist(i, 2) = data(end,3); % earth X
   fordist(i, 3) = data(end,4); % earth Y
   fordist(i, 4) = data(end,11); % apollo X
   fordist(i, 5) = data(end,12); % apollo Y
   

%    forvmax(i, 6) = data(1,12); % apollo Y initial
   
   t(1,i) = data(end, 1);
end

%% COMPUTING DISTANCE
RT = 6378.1 * 1000; % earth's radius

dist = zeros(nsimul, 1);
r0 = zeros(nsimul, 1);
for i=1:nsimul
   dist(i, 1) = sqrt((fordist(i,4) - fordist(i,2))^2 + (fordist(i, 5) - fordist(i,3))^2 ) - RT;
end

%% PLOTTING DATA
fig1=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(epsilon, dist,'x');
grid on;

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel("Precision $\epsilon$ [m]");
ylabel("Distance between Earth and Apollo 13 [m]");

hold off;