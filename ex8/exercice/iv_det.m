%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
sigma_norm = 0.06;
n = 11.06;
tfin = 5000;
Ninters = 300;
tdetect = 1000;

dt = 1.;

dL = Ninters*xL/(xL-xR);
dR = Ninters;


% Choisir le bon delta
% delta = 10 % E > v0
Etmp=0.0257;
% delta = sqrt(2*Etmp/(omega^2)); % E = V0
% delta = 150 % E < v0
delta = 64;

x0 = -delta;

nsimul = 1;


output = sprintf("iv_det");



cmd = sprintf("./Exercice8 configuration.in output=%s xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f t_detect=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f dL=%0.15f dR=%0.15f", output, xL, xR, omega, delta, x0, sigma_norm, tdetect, n, tfin, Ninters, dt, dL, dR);




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

%Evo

% figEvo=figure;
% hold on;
% 
% set(gca, 'fontsize', 25);
% set(gca, 'LineWidth',1.5);
% 
% pcolor(X, T, psi2);
% %colbar=colorbar;
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
% tick2 = changePrecision(mmin+div10, 2);
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
% %Prob
% 
% figProb=figure;
% hold on;
% 
% set(gca, 'fontsize', 25);
% set(gca, 'LineWidth',1.5);
% 
% plot(t, probG, '-', 'linewidth', 2);
% plot(t, probD, '-', 'linewidth', 2);
% 
% box on
% grid on;
% 
% xlabel("$t~[t_P]$");
% ylabel("$P(t)$");
% 
% legend(["$P_{x<0}(t)$", "$P_{x>0}(t)$"])
% 
% hold off;

%xmoy

figxmoy=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, xmoy, '-', 'linewidth', 1.2);

box on
grid on;

xlabel("$t~[t_P]$");
ylabel("$\langle x \rangle (t)~[\ell_P]$");

hold off;

%pmoy

figpmoy=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, pmoy, '-', 'linewidth', 1.2);

box on
grid on;

xlabel("$t~[t_P]$");
ylabel("$\langle p \rangle (t)~[m_P~c]$");

hold off;

%E

figH=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, E, '-', 'linewidth', 1.2);

box on
grid on;

xlabel("$t~[t_P]$");
ylabel("$\langle H \rangle (t)~[E_P]$");

hold off;

%delx

figdelx=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, delx, '-', 'linewidth', 1.2);

box on
grid on;

xlabel("$t~[t_P]$");
ylabel("$\langle \Delta x \rangle (t)~[\ell_P]$");

hold off;

%delp

figdelp=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(t, delp, '-', 'linewidth', 1.2);

box on
grid on;

xlabel("$t~[t_P]$");
ylabel("$\langle \Delta p \rangle (t)~[m_P~c]$");

hold off;

%psi2

figpsi=figure;
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(x,psi2(1001, :),'linewidth',1.2);
box on;
grid on;

xlabel("$x~[\ell_P]$");
ylabel("$|\psi(x,t=3000)|^2$");
hold off;

saveas(figxmoy,'graphs/iv_det_xmoy','epsc');
saveas(figpmoy,'graphs/iv_det_pmoy','epsc');
saveas(figH,'graphs/iv_det_H','epsc');
saveas(figdelx,'graphs/iv_det_delx','epsc');
saveas(figdelp,'graphs/iv_det_delp','epsc');

%% Function
function n = changePrecision(value, digits)
    if value ~= 0
        order = floor(log(abs(value))./log(10)) - digits + 1; % inspired by https://ch.mathworks.com/matlabcentral/fileexchange/28559-order-of-magnitude-of-number
        n = floor(value*10^-order)*10^order;
    else
        n = 0;
    end
end