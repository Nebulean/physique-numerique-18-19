%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
sigma_norm = 0.06;
n = 11;
tfin = 5000;
Ninters = 300;

dt = 1;


% Choisir le bon delta
% delta = 10 % E > v0
Etmp=0.0257;
% delta = sqrt(2*Etmp/(omega^2)); % E = V0
% delta = 150 % E < v0
delta = 64;

x0 = -delta;

nsimul = 1;


output = sprintf("iv_neg");



cmd = sprintf("./Exercice8 configuration.in output=%s xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", output, xL, xR, omega, delta, x0, sigma_norm, n, tfin, Ninters, dt);




%% Simulations

disp(cmd);
system(cmd);



%% Taitement des données
% obs: t, probD, probG, E, xmoy, x2moy, pmoy, p2moy
psi2 = load(sprintf("%s_psi2.out", output));

data = load(sprintf("%s_obs.out", output))
t = data(:,1);
probG = data(:,2);
probD = data(:,3);
E = data(:,4);
xmoy = data(:,5);
x2moy = data(:,6);
pmoy = data(:,7);
p2moy = data(:,8);
delx = data(:,9);
delp = data(:,10);

data = load(sprintf("%s_pot.out", output));
x = data(:,1);
V = data(:,2);

[X, T] = meshgrid(x,t);


%% Figures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

figEvo=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

pcolor(X, T, psi2);
colbar=colorbar;
shading interp;

box on

xlabel("$x~[\ell_P]$");
ylabel("$t~[t_P]$");
ylabel(colbar, "$|\psi(x,t)|^2$", 'interpreter', 'latex', 'fontsize', 25);

hold off;


figProb=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, probG, '-', 'linewidth', 2);
plot(t, probD, '-', 'linewidth', 2);

box on
grid on;

xlabel("$t~[t_P]$");
ylabel("$P(t)$");

legend(["$P_{x<0}(t)$", "$P_{x>0}(t)$"])

hold off;

%psi2

figpsi=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(x,psi2(3000, :),'linewidth',1.2);
box on;
grid on;

xlabel("$x~[\ell_P]$");
ylabel("$|\psi(x,t=3000)|^2$");
hold off;