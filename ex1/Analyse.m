% Nom du fichier d'output a analyser
%filename = 'output.out'
%filename = '1000_steps.out';
%filename = '2000_steps.out';
%filename = '4000_steps.out';

%filename = '200_atm.out'
%filename = '400_atm.out'
%filename = '800_atm.out'
%filename = '1600_atm.out'
%filename = '3200_atm.out'
filename = '6400_atm.out'

%fb200 = '200_atm.out'
%fb400 = '400_atm.out'
%fb800 = '800_atm.out'
%fb1600 = '1600_atm.out'
%fb3200 = '3200_atm.out'
%fb6400 = '6400_atm.out'

%filename = 'test.out'
% Chargement des donnees
data = load(filename);

% Extraction des quantites d'interet
% (Le code c++ ecrit t, z(t) et v(t) en colonnes.)
t = data(:,1);
z = data(:,2);
v = data(:,3);

% Figures
figure('NumberTitle', 'Off', 'Name', [filename ': z(t)'])
plot(t, z, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename ': v(t)'])
plot(t, v, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename ': v(z)'])
plot(z, v, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on



%% Voici un exemple pour les etudes de convergences:
% nsteps = [1000 2000 4000 8000 16000 32000];
% zfin = [...];
% figure
% plot(nsteps, zfin, '+')
% xlabel('Nombre d''iterations')
% ylabel('z(t_{fin})')
% grid on
