%% prédéclaration
f = {};
d = {};
t = {};
z = {};
v = {};

%% CHOIX: Nom des légendes
% exercice A
% StepsName = {'1000', '2000', '4000', '8000', '16000', '32000'};

%exercice B
StepsName = {'200', '400', '800', '1600', '3200', '6400'};


%% CHOIX: chargement des fichiers les fichiers à charger
% exercice A
%f{1} = '1000A.out';
%f{2} = '2000A.out';
%f{3} = '4000A.out';
%f{4} = '8000A.out';
%f{5} = '16000A.out';
%f{6} = '32000A.out';

% exercice B
f{1} = '200B.out';
f{2} = '400B.out';
f{3} = '800B.out';
f{4} = '1600B.out';
f{5} = '3200B.out';
f{6} = '6400B.out';

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
figure('NumberTitle', 'Off', 'Name', 'z(t)');
hold on; % permet de dessiner plusieurs trucs
xlabel('t [s]');
ylabel('z [m]');
grid on;

% dessine les résultats
for i = 1:length(f)
    plot(t{i}, z{i}, '-');
end

% Un peu de déco
legend(StepsName); % donne des noms aux lignes
legend('boxoff'); % enlève la boite autours des légendes
title(legend, "steps"); % ajoute un titre aux légendes

hold off; % de la figure


%% Graph de convergence
% exercice A
%nsteps = [1000 2000 4000 8000 16000 32000];
% exercice B
%nsteps = [200 400 800 1600 3200 6400];

%for i=1:length(f)
%    zfin(i) = z{i}(length(z{i}));
%end
%figure;
%plot(nsteps, zfin, '-');
%hold on;
%plot(nsteps, zfin, '+');
%hold off;
%xlabel('Nombre d''iterations');
%ylabel('z(t_{fin})');
%grid on;
