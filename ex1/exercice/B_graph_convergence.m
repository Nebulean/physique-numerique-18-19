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
nsteps = [200 400 800 1600 3200 6400];

for i=1:length(f)
    zfin(i) = z{i}(length(z{i}));
end
figure;
plot(nsteps, zfin, '-');
hold on;
plot(nsteps, zfin, '+');
hold off;
xlabel('Nombre d''itérations');
ylabel('z( t = 10s ) [m]');
grid on;

%% Graph de convergence de v
nsteps = [200 400 800 1600 3200 6400];

for i=1:length(f)
    vfin(i) = v{i}(length(v{i}));
end
figure;
plot(nsteps, vfin, '-');
hold on;
plot(nsteps, vfin, '+');
hold off;
xlabel('Nombre d''itérations');
ylabel('v( t = 10s ) [m/s]');
grid on;
