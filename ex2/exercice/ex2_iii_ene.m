%% On commence par dessiner les trajectoires
repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

output='ex2_iii_energy_ALL'
schema = ["Euler", "EulerCromer", "RungeKutta2"];

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

darkblue = [0.2 0.2 0.8];
darkred = [0.8 0.2 0.2];
darkgreen = [0.2 0.8 0.2];
darkyellow = [0.8 0.8 0.2];
% Copie direct des valeurs de configuration.in
q      = 1.6022e-19;
m      = 1.6726e-27;
B0     = 3;
E      = 0;
Kappa  = 0;
vx0    = 0;
vy0    = 4e5;
x0     = -vy0*m/(q*B0);
y0     = 0;

v0 = vy0;
omega=q*B0/m;

for i=1:3
% Execution du programme en lui envoyant la valeur a scanner en argument
cmd = sprintf('./Exercice2 configuration.in output=%s-%s.out schema=%s nsteps=1000', output, schema(i), schema(i));
disp(cmd)
system(cmd);
end

f1 = figure
hold on
t = [];
color = [darkblue; darkred; darkgreen];
for i=1:3
    filename = sprintf('%s-%s.out', output, schema(i));
    d = load(filename);
    t = d(:,1); % Le t reste constant
    energy = d(:,6);

    plot(t,energy, 'LineWidth',1, 'Color', color(i,:));
end

xlabel('$t$ [s]');
ylabel('Energy [J]');
set(gca, 'fontsize',20);

legstr = ["Euler", "Euler Cromer", "Runge Kutta 2"];

l = legend(legstr);
set(l, 'Location', 'northwest');
set(l, 'FontSize', 18);

grid on

hold off

saveas(f1, 'graphs/ex2_iii_ene', 'epsc');




% %% On commence par dessiner les trajectoires
% repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
% executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
% input = 'configuration.in'; % Nom du fichier d'entree de base
% 
% output='ex2_iii_energy_NoEuler'
% schema = ["Euler", "EulerCromer", "RungeKutta2"];
% 
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'defaultLegendInterpreter', 'latex');
% set(groot, 'defaultTextInterpreter', 'latex');
% set(groot, 'defaultAxesFontSize', 18);
% 
% darkblue = [0.2 0.2 0.8];
% darkred = [0.8 0.2 0.2];
% darkgreen = [0.2 0.8 0.2];
% darkyellow = [0.8 0.8 0.2];
% % Copie direct des valeurs de configuration.in
% q      = 1.6022e-19;
% m      = 1.6726e-27;
% B0     = 3;
% E      = 0;
% Kappa  = 0;
% vx0    = 0;
% vy0    = 4e5;
% x0     = -vy0*m/(q*B0);
% y0     = 0;
% 
% v0 = vy0;
% omega=q*B0/m;
% 
% for i=2:3
% % Execution du programme en lui envoyant la valeur a scanner en argument
% cmd = sprintf('./Exercice2 configuration.in output=%s-%s.out schema=%s nsteps=1000', output, schema(i), schema(i));
% disp(cmd)
% system(cmd);
% end
% 
% f2 = figure
% hold on
% t = [];
% color = [darkblue; darkred; darkgreen];
% for i=2:3
%     filename = sprintf('%s-%s.out', output, schema(i));
%     d = load(filename);
%     t = d(:,1); % Le t reste constant
%     energy = d(:,6);
% 
%     plot(t,energy, 'LineWidth',1, 'Color', color(i,:));
% end
% 
% xlabel('$t$ [s]');
% ylabel('Energy [J]');
% set(gca, 'fontsize',20);
% 
% legstr = ["Euler Cromer", "Runge Kutta 2"];
% 
% l = legend(legstr);
% set(l, 'Location', 'northwest');
% set(l, 'FontSize', 18);
% 
% grid on
% 
% hold off