%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
delta = 0;
x0 = 0;
sigma_norm = 0.06;
n = 14;
tfin = 5000;
Ninters = 300;

init2 = 7;
end2 = 15;
dt = tfin./2.^(init2:end2);

nsimul = length(dt);

output = {};
for i=1:nsimul
    output{i} = sprintf("i_conv_%i", i);
end

cmd = {};
for i=1:nsimul
    cmd{i} = sprintf("./Exercice8 configuration.in output=%s xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", output{i}, xL, xR, omega, delta, x0, sigma_norm, n, tfin, Ninters, dt(i));
end




%% Simulations
for i=1:nsimul
    disp(cmd{i});
    %system(cmd{i});
end





%% Traitement des données.
% osb: t, probD, probG, E, xmoy, x2moy, pmoy, p2moy
for i=1:nsimul
    data = load(sprintf("%s_obs.out", output{i}));
    
    N(i) = length(data(:,1));
    xmoy(i) = data(end, 5);
    pmoy(i) = data(end, 7);
    delx(i) = data(end, 9);
    delp(i) = data(end, 10);
end

%



%% Figures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

figx = figure;
hold on;

X = N;
Y = abs(xmoy-xmoy(end));
X(end) = [];
Y(end) = [];

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(X, Y, 'x', 'markersize', 10, 'linewidth', 1.5);

[fit, slope] = poly_approx(X, Y, 1, 2, true);

plot(fit(:,1), fit(:,2), '-', 'linewidth', 1.5);

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel("$N$");
ylabel("$|\langle x \rangle - \langle x_{best}\rangle|~[\ell_P]$");

legend(["data", sprintf("slope = %0.5f", slope)]);
grid on;
box on;


hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figp = figure;
hold on;

X = N;
Y = abs(pmoy-pmoy(end));
X(end) = [];
Y(end) = [];

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(X, Y, 'x', 'markersize', 10, 'linewidth', 1.5);

[fit, slope] = poly_approx(X, Y, 1, 2, true);

plot(fit(:,1), fit(:,2), '-', 'linewidth', 1.5);

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel("$N$");
ylabel("$|\langle p \rangle - \langle p_{best}\rangle|~[m_P~c]$");

legend(["data", sprintf("slope = %0.5f", slope)]);
grid on;
box on;


hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figdx = figure;
hold on;

X = N;
Y = abs(delx-delx(end));
X(end) = [];
Y(end) = [];

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(X, Y, 'x', 'markersize', 10, 'linewidth', 1.5);

[fit, slope] = poly_approx(X, Y, 1, 2, true);

plot(fit(:,1), fit(:,2), '-', 'linewidth', 1.5);

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel("$N$");
ylabel("$|\langle \Delta x \rangle - \langle \Delta x_{best}\rangle|~[\ell_P]$");

legend(["data", sprintf("slope = %0.5f", slope)]);
grid on;
box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figdp = figure;
hold on;

X = N;
Y = abs(delp-delp(end));
X(end) = [];
Y(end) = [];

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(X, Y, 'x', 'markersize', 10, 'linewidth', 1.5);

[fit, slope] = poly_approx(X, Y, 1, 2, true);

plot(fit(:,1), fit(:,2), '-', 'linewidth', 1.5);

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel("$N$");
ylabel("$|\langle \Delta p \rangle - \langle \Delta p_{best}\rangle|~[m_P~c]$");

legend(["data", sprintf("slope = %0.5f", slope)]);
grid on;
box on;


hold off;


%% Saves
saveas(figx, "graphs/i_conv_x", "epsc");
saveas(figp, "graphs/i_conv_p", "epsc");
saveas(figdx, "graphs/i_conv_dx", "epsc");
saveas(figdp, "graphs/i_conv_dp", "epsc");

%% Fonctions
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