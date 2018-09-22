%% prédéclaration
f = {};
d = {};
t = {};
z = {};
v = {};

%% VARIABLES A CHANGER
XAxisName = 't [s]';
YAxisName = 'z [m]';
%StepsName = {'200', '400', '800', '1600', '3200', '6400'};
StepsName = {'1000', '2000', '4000', '8000', '16000', '32000'};

%% On choisi les fichiers à charger
% exercice A
f{1} = '1000A.out';
f{2} = '2000A.out';
f{3} = '4000A.out';
f{4} = '8000A.out';
f{5} = '16000A.out';
f{6} = '32000A.out';

% exercice B
%f{1} = '200B.out';
%f{2} = '400B.out';
%f{3} = '800B.out';
%f{4} = '1600B.out';
%f{5} = '3200B.out';
%f{6} = '6400B.out';

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
figure('NumberTitle', 'Off', 'Name', 'graph');
plot(t{1}, z{1}, '-');
xlabel(XAxisName);
ylabel(YAxisName);
grid on;
hold on; % permet de dessiner d'autres trucs dans la figure

for i = 2:length(f)
    plot(t{i}, z{i}, '-');
end

% Un peu de déco
legend(StepsName); % donne des noms aux lignes
legend('boxoff'); % enlève la boite autours des légendes
title(legend, "steps"); % ajoute un titre aux légendes


hold off; % de la figure
