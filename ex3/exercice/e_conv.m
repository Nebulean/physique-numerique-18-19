%clear

%% PREDECLARATIONS
repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

% SI LES SIMULATIONS ET CALCULS SONT DEJA FAITES
simul_done = 1; % 0 = A FAIRE, 1 = DEJA FAIT
calculs_done = 1; % 0 = A FAIRE, 1 = DEJA FAIT

g = 9.81;
L = 0.1;
omega0 = sqrt(g/L)

Omega = omega0
d = 0.04
kappa = 0.
theta0 = 0.
thetadot0 = 1e-2

tfin = 20*(2*pi)/Omega
%tfin = 12.6874797;
nsimul = 20;

%nsteps = 2.^(14:18);
nsteps = round(logspace(4,5,nsimul));

dt = tfin ./ nsteps;

%% SIMULATIONS
output = cell(1, nsimul);
for i=1:nsimul
    output{i} = sprintf('dt=%f.out', dt(i));
    
    cmd = sprintf('%s%s %s Omega=%0.15f d=%0.15f kappa=%0.15f theta0=%0.15f thetadot0=%0.15f tFin=%0.15f dt=%0.15f sampling=1 output=%s', repertoire, executable, input, Omega, d, kappa, theta0, thetadot0, tfin, dt(i), output{i});
    
    disp(cmd);
    if simul_done == 0
        system(cmd);
    end
    
    fprintf('Simulation actuelle: %d\n', i);
end

%% CALCULS
if calculs_done == 0
    thetafin = zeros(1,nsimul);
    tfinal = zeros(1,nsimul);
    for i=1:nsimul
        d = load(output{i});
        
        tfinal(i) = d(end,1)
        thetafin(i) = d(end, 2)
    end
    dt2 = dt.^2; % ça marche pas avec dt^2, j'comprend pas :(
end

%% PLOT
f=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(dt2, thetafin, 'x');

ylim([-6.988e-3, -6.968e-3]);

xlabel("$\Delta t^2$");
ylabel("$\theta\left(t_{end}\right)$");

polynome = poly_approx(dt2(1,1:5), thetafin(1,1:5), 1, 500);

fit = polyfit(dt2(1,1:5), thetafin(1,1:5), 1);
sprintf('%0.15f', fit(1))

plot(polynome(1,:), polynome(2,:),'-');
grid on;

hold off;

saveas(f, 'graphs/e_conv_dt2','epsc');




function polynome = poly_approx(x, y, ordre, steps)
    pf = polyfit(x, y, ordre);
    T = linspace(min(x), max(x), steps);
    
    n = ordre + 1;
    
    polynome = zeros(2,length(T));
    for i=1:n
       polynome(2,:) = polynome(2,:) + pf(i)*T.^(n-i);
    end
    polynome(1,:) = T;
end