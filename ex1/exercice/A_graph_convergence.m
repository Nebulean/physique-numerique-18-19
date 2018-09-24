%% prédéclaration
f = {};
d = {};
t = {};
z = {};
v = {};

%% Enregistrement du nom des fichiers
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

%% Graph de convergence de z
nstepsZ = [1000 2000 4000 8000 16000 32000];

for i=1:length(f)
    zfin(i) = z{i}(length(z{i}));
end
figure;
plot(nstepsZ, zfin, '-');
hold on;
plot(nstepsZ, zfin, '+');
hold off;
xlabel('Nombre d''itérations');
ylabel('z( t = 86400s ) [m]');
grid on;

%% Graph de convergence de v
nstepsV = [1000 2000 4000 8000 16000 32000];

for i=1:length(f)
    vfin(i) = v{i}(length(v{i}));
end
figure;
plot(nstepsV, vfin, '-');
hold on;
plot(nstepsV, vfin, '+');
hold off;
xlabel('Nombre d''itérations');
ylabel('v( t = 86400s ) [m/s]');
grid on;
