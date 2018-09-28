%% prédéclaration
f = {};
d = {};
t = {};
z = {};
v = {};

%% Nom des légendes
StepsName = {'1000', '2000', '4000', '8000', '16000', '32000'};


%% chargement des fichiers les fichiers à charger
f{1} = '1000A.out';
f{2} = '2000A.out';
f{3} = '4000A.out';
f{4} = '8000A.out';
f{5} = '16000A.out';
f{6} = '32000A.out';

%% On traite les données
for i = 1:length(f)
    %% on commence par charger les données
    d{i} = load(f{i});

    %% ensuite on sépare les colonnes
    t{i} = d{i}(:,1);
    z{i} = d{i}(:,2);
    v{i} = d{i}(:,3);
end

%% Maintenant on va créer les graphiques
% Graphique de vitesse
figure('NumberTitle', 'Off', 'Name', 'v(t)');
hold on; % Début de la figure
xlabel('t [s]');
ylabel('v [m]');
set(gca,'fontsize',16); % Pour changer la taille de la police des axes
set(gca, 'LineWidth', 1.2); % Rend les axes plus gras
grid on;

% dessine les résultats
for i = 1:length(f)
    p = plot(t{i}, v{i}, '-');
    p.LineWidth = 1.5; % Pour rendre les lignes plus grasses.
end

% Ajout d'une légande
legend(StepsName); % donne des noms aux lignes
legend('boxoff'); % enlève la boite autours des légendes
title(legend, "steps"); % ajoute un titre aux légendes
set(legend, 'FontSize', 16); % Pour changer la taille de la police de la légende

hold off; % fin de la figure

% Graphique de position
figure('NumberTitle', 'Off', 'Name', 'z(t)');
hold on; % Début de la figure
xlabel('t [s]');
ylabel('z [m]');
set(gca,'fontsize',16); % Pour changer la taille de la police des axes
set(gca, 'LineWidth', 1.2); % Rend les axes plus gras
grid on;

% dessine les résultats
for i = 1:length(f)
    p = plot(t{i}, z{i}, '-');
    p.LineWidth = 1.5; % Pour rendre les lignes plus grasses.
end

% Ajout d'une légande
legend(StepsName); % donne des noms aux lignes
legend('boxoff'); % enlève la boite autours des légendes
title(legend, "steps"); % ajoute un titre aux légendes
set(legend, 'FontSize', 16); % Pour changer la taille de la police de la légende

hold off; % fin de la figure