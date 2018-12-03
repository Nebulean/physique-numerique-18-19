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
fordist = zeros(nsimul, 7);
tfin = zeros(1, nsimul);
forvmax = zeros(nsimul, 6);
for i=1:nsimul
   data = load(output{1,i});
   
   fordist(i, 1) = data(end,2); % dt
   fordist(i, 2) = data(end,3); % earth X
   fordist(i, 3) = data(end,4); % earth Y
   fordist(i, 4) = data(end,11); % apollo X
   fordist(i, 5) = data(end,12); % apollo Y
   
   forvmax(i, 1) = data(1,13); % apollo vx initial
   forvmax(i, 2) = data(1,14) % apollo vy initial
   forvmax(i, 3) = data(1,3); % earth X initial
   forvmax(i, 4) = data(1,4); % earth Y initial
   forvmax(i, 5) = data(1,11); % apollo X initial
   forvmax(i, 6) = data(1,12) % apollo Y initial
   
   t(1,i) = data(end, 1);
end

%% COMPUTING DISTANCE
RT = 6378.1 * 1000; % earth's radius

dist = zeros(nsimul, 1);
vmax = zeros(nsimul, 1);
r0 = zeros(nsimul, 1);
for i=1:nsimul
   dist(i, 1) = sqrt((fordist(i,4) - fordist(i,2))^2 + (fordist(i, 5) - fordist(i,3))^2 ) - RT;
   vmax(i, 1) = sqrt(fordist(i,6)^2 + fordist(i,7)^2);
   r0(i, 1) = sqrt( (forvmax(i,4) - forvmax(i,2))^2 + (forvmax(i,5) - forvmax(i,3))^2 );
end


G = 6.67408e-11; % gran constant
M = 5.972e24; % mass of earth

vmaxTH = sqrt(( forvmax(1, 1)^2 + forvmax(1,2)^2 ) + G*M/2 .* abs(1./dist - 1./r0) );


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