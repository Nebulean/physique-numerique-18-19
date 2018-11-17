%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 101; % Nombre de simulations a faire

%Copie des paramètres de configuration.in
g = 9.81;
L = 0.1;
dt = 0.02;

% Extended plot
omega0 = sqrt(g/L);
epsilon = 9*10^(-1);
Omega = linspace(2*omega0-epsilon,2*omega0+epsilon,nsimul);

% close to explosion
% omega0 = 19.15;
% epsilon = 0.1*10^(-1);
% Omega = linspace(omega0-epsilon,omega0+epsilon,nsimul);

paramstr = 'Omega';
param = Omega;

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i),6), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s dt=%s theta0=0. thetadot0=1e-2 d=0.005 kappa=0.05 tFin=100', repertoire, executable, input, paramstr, param(i), output{i}, dt);
    disp(cmd)
    %system(cmd);
end


%% Analyse %%
%%%%%%%%%%%%%

Emax=zeros(1,nsimul);

for i=1:nsimul
    data = load(output{i});
    emec = data(:,4);
    Emax(i)=max(emec);
end

%% Figures %%
%%%%%%%%%%%%%

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

f1=figure;
hold on
grid on
plot(Omega,Emax,'x')

% Pour le close plot
% xlim([19.139, 19.161]);

set(gca, 'fontsize', 22);
xlabel('$\Omega [rad/s]$')
ylabel('$\max_t E_{mec}(t) [J]$')

%saveas(f1, 'graphs/deletethis', 'epsc') % -> sauvegarde manuelle vu que
%y'a deux graphs.