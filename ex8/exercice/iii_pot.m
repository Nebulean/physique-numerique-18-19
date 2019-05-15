%% Paramètres et initialisation
xL = -200;
xR = 200;
omega = 0.003;
sigma_norm = 0.06;
n = 14;
tfin = 5000;
Ninters = 300;

dt = 1;

Etmp=0.0257; % Marche que dans un cas. Changer si les paramètres changement.
% Choisir le bon delta
% delta = 10 % E > v0
% delta = sqrt(2*Etmp/(omega^2)); % E = V0
delta = [10, sqrt(2*Etmp/(omega^2)), 150] % E < v0
%delta = 100;

% n = 30% test

x0 = -delta;

nsimul = 3;

output = {};
for i=1:nsimul
    output{i} = sprintf("iii_evo_%i", i);
end

cmd = {};
for i=1:nsimul
    cmd{i} = sprintf("./Exercice8 configuration.in output=%s xL=%0.15f xR=%0.15f omega=%0.15f delta=%0.15f x0=%0.15f sigma_norm=%0.15f n=%0.15f tfin=%0.15f Ninters=%0.15f dt=%0.15f", output{i}, xL, xR, omega, delta(i), x0(i), sigma_norm, n, tfin, Ninters, dt);
end



%% Simulations
for i=1:nsimul
    disp(cmd{i});
    system(cmd{i});
end

%% Traitement des données
V = {};
for i=1:nsimul
    data = load(sprintf("%s_pot.out", output{i}));
    x = data(:,1);
    V{i} = data(:,2);
    
    data = load(sprintf("%s_obs.out", output{i}))
    E(i) = mean(data(:,4));
end


%% Figure
figpot=figure; % Pas dans le rapport
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

for i=1:nsimul
    plot(x, V{i}, '-', 'linewidth', 2);
    % line([xL, xR], [E(i) E(i)], 'linewidth', 2);
end


xlabel("$x~[\ell_P]$");
ylabel("$V~[E_P]$");

legend("$\Delta = 10$", "$\Delta = 75.57$", "$\Delta = 150$");

box on;
grid on;

hold off;