%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
sigma_norm = 0.06;
n = 14;
tfin = 10000;
Ninters = 2500;
delta = 64;
dt = 5;
x0 = -delta;

nsimul = 21;

% delta = logspace(-2, 1, nsimul);
delta = linspace(0, 6, nsimul);

output = {};
for i=1:nsimul
   output{i} = sprintf("v_bonus_conv_pot_delta=%0.5f", delta(i));
end

%% Simulations
for i=1:nsimul
    cmd = sprintf("./Exercice8 configuration.in output=%s xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", output{i}, xL, xR, omega, delta(i), x0, sigma_norm, n, tfin, Ninters, dt);
    disp(cmd);
	% system(cmd);
    cmd = "";
end

%% Traitement des données
% obs: t, probD, probG, E, xmoy, x2moy, pmoy, p2moy
y_axis1 = zeros(nsimul,1);
y_axis2 = zeros(nsimul,1);
E0 = zeros(nsimul,1);
for i=1:nsimul
    data = load(sprintf("%s_obs.out", output{i}));
    t = data(:,1);
    probG = data(:,2);
    probD = data(:,3);
    E = data(:,4);
    E0(i) = E(1);
    
    [time, Tidx] = min(abs(t-1000));
    
    y_axis1(i,1) = probG(Tidx);
    y_axis2(i,1) = probD(Tidx);
end

psi2 = load(sprintf("%s_psi2.out", output{end}));
% 
% data = load(sprintf("%s_obs.out", output))
% t = data(:,1);
% probG = data(:,2);
% probD = data(:,3);
% E = data(:,4);
% xmoy = data(:,5);
% x2moy = data(:,6);
% pmoy = data(:,7);
% p2moy = data(:,8);
% delx = data(:,9);
% delp = data(:,10);
% 
data = load(sprintf("%s_pot.out", output{end}));
x = data(:,1);
V = data(:,2);
% 
[X, T] = meshgrid(x,t);


%% Figures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

figConv=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

% plot(delta, y_axis1, 'x', 'markersize', 10, 'linewidth', 1.5);
plot(delta, y_axis2, 'x', 'markersize', 10, 'linewidth', 1.5, 'markersize', 10);
% approx = poly_approx(delta(3:11)', y_axis2(3:11), 1, 1000, false);
% plot(approx(:,1), approx(:,2), '-');

[expfit, pf] = poly_approx(delta(1:end)', log(y_axis2(1:end)), 1, 1000, false);
plot(expfit(:,1), exp(expfit(:,2)), '--', 'linewidth', 1.5);


% set(gca, 'xscale', 'log');
% set(gca, 'yscale', 'log');
% set(gca,'yTickLabel',{'e^-1.28','e^-1','e^0','e^1'})

legend(["data", sprintf("$y=%0.3fe^{%0.3fx}$", exp(pf(2)), pf(1))])


xlabel("$\Delta~[\ell_P]$");
ylabel("$P_{x>0}$")

ylim([0, 1]);

grid on;
box on;

hold off;


% test=figure
% hold on;
% m=1;
% hbar=1;
% % E0 = E(1);
% V0 = 0.1;
% alpha = sqrt(2*m*(V0-E(1)))/hbar;
% a = 2.*delta;
% 
% T = (1+a'.*sinh(alpha)^2./(4.*E0./V0.*(1-E0./V0))).^(-1)
% 
% plot(delta, T, 'x');
% plot(delta, y_axis2, 'x')
% 
% hold off;

% figEvo=figure;
% hold on;
% 
% set(gca, 'fontsize', 25);
% set(gca, 'LineWidth',1.5);
% 
% pcolor(X, T, psi2);
% 
% colbar=colorbar('TickLabelInterpreter', 'latex', 'fontsize', 25);
% colbar.Label.String = '$|\psi(x,t)|^2$';
% colbar.Label.Interpreter = 'latex';
% mmax = max(max(psi2));
% mmin = min(min(psi2));
% div10 = (mmax-mmin)/4;
% div2 = (mmax-mmin)/2;
% colbar.Label.Position = [3, mmin+div2, 0];
% tick1 = changePrecision(mmin, 2);
% tick2 = changePrecision(mmin+div10, 1);
% tick3 = changePrecision(mmax-div10, 2);
% tick4 = changePrecision(mmax, 2);
% set(colbar, 'ytick', [tick1, tick2, tick3, tick4]);
% 
% shading interp;
% 
% box on
% 
% xlabel("$x~[\ell_P]$");
% ylabel("$t~[t_P]$");
% ylabel(colbar, "$|\psi(x,t)|^2$", 'interpreter', 'latex', 'fontsize', 25);
% 
% hold off;


% 
% figProb=figure;
% hold on;
% 
% set(gca, 'fontsize', 25);
% set(gca, 'LineWidth',1.5);
% 
% plot(t, probG, '-', 'linewidth', 1.5);
% plot(t, probD, '-', 'linewidth', 1.5);
% 
% xlabel("$t~[t_P]$");
% ylabel("P(t)");
% 
% box on;
% grid on;
% 
% legend("$P_{x<0}(t)$", "$P_{x>0}(t)$");
% 
% hold off;
% 

figpot=figure; % Pas dans le rapport
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(x, V, '-', 'linewidth', 1.5, 'color', 'blue');
line([xL, xR], [E(1) E(1)], 'color', 'red', 'linewidth', 1.5);

xlabel("$x~[\ell_P]$");
ylabel("$V~[E_P]$");

legend(["$V(x)$", "E"], 'location', 'northeast');

box on;
grid on;

hold off;


%% Save
saveas(figConv, "graphs/v_conv_delta", "epsc");
saveas(figpot, "graphs/v_conv_delta_pot", "epsc");

%% Function
function n = changePrecision(value, digits)
    if value ~= 0
        order = floor(log(abs(value))./log(10)) - digits + 1; % inspired by https://ch.mathworks.com/matlabcentral/fileexchange/28559-order-of-magnitude-of-number
        n = floor(value*10^-order)*10^order;
    else
        n = 0;
    end
end


function [polynome, param] = poly_approx(x, y, ordre, steps, isLog)
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
    length(x)
    length(y)
    pf = polyfit(x, y, ordre);
    %coefvalues(pf)
    param = pf;
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

