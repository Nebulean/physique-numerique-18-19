%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 5; % Nombre de simulations a faire

%Copie des paramètres de configuration.in
g = 9.81;
L = 0.1;
dt = 0.02;

omega0 = sqrt(g/L);
epsilon = 10^(-1); %variation de Omega

Omega = linspace(omega0-epsilon,omega0+epsilon,nsimul);

paramstr = 'Omega';
param = Omega;

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s dt=%s theta0=0. thetadot0=1e-2 d=0.03 kappa=0. tFin=250', repertoire, executable, input, paramstr, param(i), output{i}, dt);
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

%Copie des paramètres de configuration.in
tfin = 250;
dt = 0.02;

thetavar=cell(1,nsimul);
thetadotvar=cell(1,nsimul);
emecvar=cell(1,nsimul);

c1 = load(output{1});
t = c1(:,1); %load le temps de la première simul (c'est le même pour toutes)

for i = 1:nsimul
    data = load(output{i});
    
    thetavar{i} = data(:,2);
    thetadotvar{i} = data(:,3);
    emecvar{i} = data(:,4);
    
end

%% Figures %%
%%%%%%%%%%%%%

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

leg = strings(1,nsimul); %vecteur de strings pour la légende

f1=figure;
hold on;
for i=1:nsimul
    plot(t, thetavar{i});
    leg(1,i) = paramstr + "=" + num2str(param(i));
end
set(gca, 'fontsize', 22);
grid on;
xlabel('$t [s]$')
ylabel('$\theta [rad]$')
legend(leg);
hold off;

f2=figure;
hold on
for i=1:nsimul
    plot(t,thetadotvar{i});
end
set(gca, 'fontsize', 22);
grid on;
xlabel('$t [s]$')
ylabel('$\dot{\theta} [rad/s]$')
legend(leg);
hold off

f3=figure;
hold on
for i=1:nsimul
    plot(t,emecvar{i});
end
set(gca, 'fontsize', 22);
grid on;
xlabel('$t [s]$')
ylabel('$E_{mec} [J]$')
legend(leg);
hold off

saveas(f1, 'graphs/c_theta', 'epsc')
saveas(f2, 'graphs/c_thetadot', 'epsc')
saveas(f3, 'graphs/c_emec', 'epsc')