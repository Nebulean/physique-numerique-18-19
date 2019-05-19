%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
sigma_norm = 0.06;
tfin = 5000;
Ninters = 300;

dt = 5;

delta = 64
x0 = -delta;

n=1:25; % Meilleur: n=11, n=1:15 pour l'évolution

nsimul = length(n)

output = {};
for i=1:nsimul
    output{i} = sprintf("iii_evo_n=%0.2f", n(i));
end

cmd = {}
for i=1:nsimul
    cmd{i} = sprintf("./Exercice8 configuration.in output=%s xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", output{i}, xL, xR, omega, delta, x0, sigma_norm, n(i), tfin, Ninters, dt);
end



%% Simulations

for i=1:nsimul
    disp(cmd{i});
    %system(cmd{i});
end


%% Taitement des données
% obs: t, probG, probD, E, xmoy, x2moy, pmoy, p2moy
for i=1:nsimul
    data = load(sprintf("%s_obs.out", output{i}));
    
    t = data(:,1);
    probG = data(:,2);
    probD = data(:,3);
    
    % On trouve la position après le premier passage.
    [time, Tidx] = min(abs(t-850));
    
    % On calcul la distance entre les deux proba
    diffProb(i) = abs(probG(Tidx) - probD(Tidx));
end

data = load(sprintf("%s_pot.out", output{11}));
x = data(:,1);

psi2 = load(sprintf("%s_psi2.out", output{11}));

data = load(sprintf("%s_obs.out", output{11}));
t = data(:,1);
probG = data(:,2);
probD = data(:,3);

[X, T] = meshgrid(x,t);


%% Figures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

figFindN = figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(n, diffProb, 'x', 'linewidth', 1.5, 'markersize', 10);

xlabel("$n$");
ylabel("$|P_{x<0} - P_{x>0}|$");

grid on;
box on;

hold off;

%%%%%%%%%

figProb=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, probG, '-', 'linewidth', 1.5);
plot(t, probD, '-', 'linewidth', 1.5);

xlabel("$t~[t_P]$");
ylabel("P(t)");

box on;
grid on;

legend("$P_{x<0}(t)$", "$P_{x>0}(t)$");

hold off;

%%%%%%%%%

figEvo=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

pcolor(X, T, psi2);

colbar=colorbar('TickLabelInterpreter', 'latex', 'fontsize', 25);
colbar.Label.String = '$|\psi(x,t)|^2$';
colbar.Label.Interpreter = 'latex';
mmax = max(max(psi2));
mmin = min(min(psi2));
div10 = (mmax-mmin)/4;
div2 = (mmax-mmin)/2;
colbar.Label.Position = [3, mmin+div2, 0];
tick1 = changePrecision(mmin, 2);
tick2 = changePrecision(mmin+div10, 1);
tick3 = changePrecision(mmax-div10, 2);
tick4 = changePrecision(mmax, 2);
set(colbar, 'ytick', [tick1, tick2, tick3, tick4]);

shading interp;

box on

xlabel("$x~[\ell_P]$");
ylabel("$t~[t_P]$");
ylabel(colbar, "$|\psi(x,t)|^2$", 'interpreter', 'latex', 'fontsize', 25);

hold off;

%% saves
saveas(figFindN, "graphs/iii_findn_n", "epsc");
saveas(figEvo, "graphs/iii_findn_evo", "epsc");
saveas(figProb, "graphs/iii_findn_prob", "epsc");

%% Function
function n = changePrecision(value, digits)
    if value ~= 0
        order = floor(log(abs(value))./log(10)) - digits + 1; % inspired by https://ch.mathworks.com/matlabcentral/fileexchange/28559-order-of-magnitude-of-number
        n = floor(value*10^-order)*10^order;
    else
        n = 0;
    end
end