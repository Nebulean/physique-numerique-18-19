%% Préparation de(s) simulation(s)
schema = "B";
omega = 2*pi/900;
n_stride = 100;
Npoints = 1000;
tfin = 10000;

param = sprintf("schema=%s n_stride=%i Npoints=%i tfin=%i omega=%f", schema, n_stride, Npoints, tfin, omega);
output = "tsunami_a";

%% Simulations
cmd = sprintf("./Exercice7 configuration_tsunami.in output=%s %s", output, param);
disp(cmd);
system(cmd);

%% Chargement des données
dataf = load(sprintf("%s_f.out", output));
datau = load(sprintf("%s_u.out", output));
dataE = load(sprintf("%s_E.out", output));

%% Traitement des données
f = dataf(:,2:end);
t = dataf(:,1);

x = datau(:,1);

% On veut suivre une unique vague, du coup, je dois trouver le maximum de
% la première vague.
[pks, indexes] = findpeaks(f(end,:));


%% On fait des graphes
fig=figure
hold on;

plot(x, f(end,:))
plot(x(indexes), f(end, indexes),'o')

hold off;