%% prédéclaration
f = {};
d = {};
t = {};
z = {};
v = {};

%% Enregistrement du nom des fichiers
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

%% Graph de convergence de z
nstepsZ = [200 400 800 1600 3200 6400];

% On veut delta t et non le nb d'étapes
for i=1:length(nstepsZ)
    nstepsZ(1,i) = duration/nstepsZ(1,i);
end

for i=1:length(f)
    zfin(i) = z{i}(length(z{i}));
end
f1 = figure;
p = plot(nstepsZ, zfin, '-');
p.LineWidth = 1.5; % Pour rendre les lignes plus grasses.
hold on;
p = plot(nstepsZ, zfin, '+');
p.LineWidth = 1.5; % Pour rendre les lignes plus grasses.
set(gca,'fontsize',16); % Pour changer la taille de la police des axes
set(gca, 'LineWidth', 1.2); % Rend les axes plus gras
hold off;
xlabel('\Delta t [s]');
ylabel('z( t = 10s ) [m]');
grid on;

%% Graph de convergence de v
nstepsV = [1000 2000 4000 8000 16000 32000];

% On veut delta t et non le nb d'étapes
for i=1:length(nstepsV)
    nstepsV(1,i) = duration/nstepsV(1,i);
end

for i=1:length(f)
    vfin(i) = v{i}(length(v{i}));
end
f2 = figure;
p = plot(nstepsV, vfin, '-');
p.LineWidth = 1.5; % Pour rendre les lignes plus grasses.
hold on;
p = plot(nstepsV, vfin, '+');
p.LineWidth = 1.5; % Pour rendre les lignes plus grasses.
set(gca,'fontsize',16); % Pour changer la taille de la police des axes
set(gca, 'LineWidth', 1.2); % Rend les axes plus gras
hold off;
xlabel('\Delta t [s]');
ylabel('v( t = 10s ) [m/s]');
grid on;

saveas(f1, "graphs/zConvB.eps", "epsc"); % exporte f1, avec le nom, en couleur (epsc).
saveas(f2, "graphs/vConvB.eps", "epsc"); % pareil