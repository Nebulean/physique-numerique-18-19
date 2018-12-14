%% TASK
% Résoudre numériquement la trajectoire d’Apollo 13 avec un pas de temps ∆t fixe. [N.B. :
% Il suffit de simuler 2 jours]. Faire une étude de convergence de l’altitude minimale hmin
% et de vmax avec ∆t et comparer avec le résultat analytique.

%% BEFORE-SIMULATIONS
tFin = 172800; % 2 days in seconds.
dtad='false'; % We do not want an adaptative dt.
atm='false';

nsimul = 10; % Number of simulations that we want. 8, because 172800 = 2^8 * 3^3 * 5^2.
% dt = zeros(nsimul, 1);
% for i=1:nsimul
%     dt(i,1) = 128/2^i; %256, 128, 64
%     fprintf("i = %d", i);
% end
dt = round(logspace(0.5, 2.5, nsimul));
%dt = round(linspace(0.1,100,nsimul));


% dt = logspace(4, 0, nsimul); % logspace ne marche pas parce que les tfin
% sont différents. Il nous faut un diviseur commun.

%% SIMULATIONS
output = cell(1, nsimul);
for i=1:nsimul
     output{1, i} = sprintf("dt=%f.out", dt(i));
     cmd = sprintf("./Exercice4 configuration.in dt=%0.15f dtad=%s tFin=%d atm=%s output=%s", dt(i), dtad, tFin, atm, output{1,i});

     disp(cmd);
     system(cmd);
end

%% TREATMENT
%% LOADING DATA
%fordist = zeros(nsimul, 7);
tfin = zeros(1, nsimul);
%val = zeros(nsimul, 6);
RT = 6378.1 * 1000; % earth's radius
distTH = 10e3 + RT;
mindist = zeros(nsimul, 1);
maxvit = zeros(nsimul, 1);
vitTH = zeros(nsimul,1);
nsteps = zeros(nsimul, 1);
h0 = zeros(nsimul, 1);
G = 6.674e-11;
M = 5.972e24;
for i=1:nsimul
   data = load(output{1,i});

%    val(i, 1) = data(end,2); % dt
%    val(i, 2) = data(end,3); % earth X
%    val(i, 3) = data(end,4); % earth Y
%    val(i, 4) = data(end,11); % apollo X
%    val(i, 5) = data(end,12); % apollo Y
%    val(i, 6) = data(end,13); % apollo vx initial
%    val(i, 7) = data(end,14) % apollo vy initial

    Ex = data(:,3);
    Ey = data(:,4);
    Ax = data(:,11);
    Ay = data(:,12);

    Avx = data(:,13);
    Avy = data(:,14);

    t = data(:, 1);
    tfin(i) = t(end);

    nsteps(i) = length(t);
    % On calcul la distance minimale
    % Distance numérique avec les dt
    dist = sqrt((Ex - Ax).^2 + (Ey - Ay).^2);
    % Premièrement, on calcul le point où il y a la distance minimale.
    [tmp, index] = min(dist);
%     if index+1 > length(dist)
%         fit = polyfit(t(index-2:index), dist(index-2:index), 2);
%     else
%         fit = polyfit(t(index-1:index+1), dist(index-1:index+1), 2);
%     end
%     A = fit(1); B = fit(2); C = fit(3);
    if index+2 > length(dist)
        fit = polyfit(t(index-4:index), dist(index-4:index), 2);
    else
        fit = polyfit(t(index-2:index+2), dist(index-2:index+2), 2);
    end
    A = fit(1); B = fit(2); C = fit(3);

    % On calcul le minimum à l'aide de l'analyse
    mindist(i, 1) = abs(C - B^2/(4*A) - distTH)% + 76.462738037109375); %TODO: Retirer les abs lorsque le résultat sera correcte
    h0(i,1) = C - B^2/(4*A) - RT;
    %% On calcul la vitesse maximale
    r0 = sqrt((Ex(1) - Ax(1))^2 + (Ey(1) - Ay(1))^2);
    % vitesse numérique avec les dt
    vit = sqrt(Avx.^2 + Avy.^2);
    vitTH(i,1) = sqrt(1200^2 + 2*G*M*(1/(RT + h0(i,1)) - 1/r0)); % Vitesse max théorique pour un hmin donné.
    [tmp, index] = max(vit);
%     if index+2 > length(vit)
%         fit = polyfit(t(index-4:index), vit(index-4:index), 2);
%     else
%         fit = polyfit(t(index-2:index+2), vit(index-2:index+2), 2);
%     end
    fit = polyfit(t(index-1:index+1), vit(index-1:index+1), 2);
    A = fit(1); B = fit(2); C = fit(3);
    maxvit(i, 1) = C - B^2/(4*A) - vitTH(i,1);
end

%% PLOTTING DATA
% hmin
fig1=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

% plot(dt, dist, 'x');
plot(dt, mindist, 'x');

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

dtfixed(:,1) = dt(1,:);
[fit, slope] = poly_approx(log10(dtfixed(1:end)), log10(mindist(1:end)), 1, 2);
fit = 10.^fit;
fitplot = plot(fit(1,:), fit(2,:), '-');

grid on;
box on;

txtReg = sprintf("slope: %0.2f", slope);
legend([fitplot], txtReg, 'Location', 'NorthWest');

xlabel("$\Delta t$ [s]");
ylabel("$|$error on $h_{min}|$ [m]");

hold off;
saveas(fig1, 'graphs/ex1b_conv_h','epsc');




% vmax
fig2=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(dt, maxvit, 'x');

grid on;
box on;

dtfixed(:,1) = dt(1,:);
[fit, slope] = poly_approx(log10(dtfixed(1:end)), log10(maxvit(1:end)), 1, 2);
fit = 10.^fit;
fitplot = plot(fit(1,:), fit(2,:), '-');

txtReg = sprintf("slope: %0.2f", slope);
legend([fitplot], txtReg, 'Location', 'NorthWest');

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel("$\Delta t$ [s]");
ylabel("$|$error on $v_{max}|$ [m/s]");

hold off;
saveas(fig2, 'graphs/ex1b_conv_vel','epsc');


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
