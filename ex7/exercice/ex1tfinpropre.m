%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

output='outpute';

n=4.
u=6.
L=20.
tfin=70.
omega=u*n*pi/L;

cmd = sprintf('%s%s %s omega=%.15g cb_droit=fixe tfin=%g output=%s', repertoire, executable, input, omega, tfin, output);
    disp(cmd);
    system(cmd);

fichier = output;
data = load([fichier,'_u.out']);
x = data(:,1);
% u = data(:,2);
data = load([fichier,'_E.out']);
t = data(:,1);
% E = data(:,2);
data = load([fichier,'_f.out']);
f = data(:,2:end);

% fth=tfin*u/L*sin(-omega/u.*x+omega*t(end))
% fth=-tfin*u/L.*cos(omega*t(end)).*sin(omega/u.*x)
fth=-tfin*u/L.*sin(omega/u.*x);

figue=figure;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

hold on;
plot(x,f(end,:));
plot(x,fth);
grid on;
box on
legend('Numerical solution','Analytical eigenmode');
xlabel('$x$ [m]')
ylabel('$f(x,t_{fin})$ [m]')
hold off;

saveas(figue,'graphs/eigenmode','epsc')