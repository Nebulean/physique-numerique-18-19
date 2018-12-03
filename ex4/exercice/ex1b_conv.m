%% TASK
% Résoudre numériquement la trajectoire d’Apollo 13 avec un pas de temps ∆t fixe. [N.B. :
% Il suffit de simuler 2 jours]. Faire une étude de convergence de l’altitude minimale hmin
% et de vmax avec ∆t et comparer avec le résultat analytique.

%% BEFORE-SIMULATIONS
tFin = 172800; % 2 days in seconds.
dtad='false'; % We do not want an adaptative dt.

nsimul = 10; % Number of simulations that we want. 8, because 172800 = 2^8 * 3^3 * 5^2.
dt = zeros(nsimul, 1);
for i=1:nsimul
    dt(i,1) = 256/2^i;
    disp(sprintf("i = %d", i));
end

% dt = logspace(4, 0, nsimul); % logspace ne marche pas parce que les tfin
% sont différents. Il nous faut un diviseur commun.

%% SIMULATIONS
output = cell(1, nsimul)
for i=1:nsimul
     output{1, i} = sprintf("dt=%f.out", dt(i));
     cmd = sprintf("./Exercice4 configuration.in dt=%0.15f dtat=%s tFin=%d output=%s", dt(i), dtad, tFin, output{1,i});
     
     disp(cmd);
     %system(cmd);
end

%% TREATMENT
%% LOADING DATA
res = zeros(nsimul, 7);
tfin = zeros(1, nsimul);
for i=1:nsimul
   data = load(output{1,i});
   
   res(i, 1) = data(end,2); % dt
   res(i, 2) = data(end,3); % earth X
   res(i, 3) = data(end,4); % earth Y
   res(i, 4) = data(end,11); % apollo X
   res(i, 5) = data(end,12); % apollo Y
   res(i, 6) = data(end,13); % apollo vx
   res(i, 7) = data(end,14) % apollo vy
   
   t(1,i) = data(end, 1);
end

%% COMPUTING DISTANCE
RT = 6378.1 * 1000; % earth's radius

dist = zeros(nsimul, 1);
vmax = zeros(nsimul, 1);
for i=1:nsimul
   dist(i, 1) = sqrt((res(i,4) - res(i,2))^2 + (res(i, 5) - res(i,3))^2 ) - RT;
   vmax(i, 1) = sqrt(res(i,6)^2 + res(i,7)^2)
end

%% PLOTTING DATA
fig1=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(dt, dist, 'x');
grid on;

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel("$\Delta t$ [s]");
ylabel("Distance between Earth and Apollo 13 [m]");

hold off;

fig2=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(dt, vmax, 'x');
grid on;

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel("$\Delta t$ [s]");
ylabel("vmax of Apollo 13 [m/s]");

hold off;