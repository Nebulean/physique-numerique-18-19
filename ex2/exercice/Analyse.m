% Nom du fichier d'output a analyser
%f = ['Euler.out', 'EulerCromer.out', 'RungeKutta2.out'];
f = {};
f{1} = 'Euler.out';
f{2} = 'EulerCromer.out';
f{3} = 'RungeKutta2.out';
f{4} = 'App2pos.out';
f{5} = 'App2neg.out';

fig = {};

output = {};
output{1} = 'graphs/graphEuler';
output{2} = 'graphs/graphEulerCromer';
output{3} = 'graphs/graphRungeKutta2';
output{4} = 'graphs/graphApp2pos';
output{5} = 'graphs/graphApp2neg';

for i=1:5
filename = f{i};
    
% Chargement des donnees
output = load(filename);

% Extraction des quantites d'interet
t = output(:,1);
x = output(:,2);
y = output(:,3);
vx = output(:,4);
vy = output(:,5);
energy = output(:,6);
mu = output(:,7);

clear output

% Figures

fig{i} = figure
subplot(2,3,1)
plot(x,y)
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')

subplot(2,3,2)
plot(vx,vy)
axis equal
grid on
xlabel('v_x [m/s]')
ylabel('v_y [m/s]')

subplot(2,3,3)
plot(t,x,t,y)
grid on
xlabel('t [s]')
ylabel('x,y [m]')
legend('x','y')

subplot(2,3,4)
p1 = plot(t,vx,t,vy)
grid on
set(p1,'LineWidth',1.5);
xlabel('t [s]')
ylabel('v_x,v_y [m/s]')
legend('v_x','Proton in constant magnetic field $B_0$')

subplot(2,3,5)
plot(t,energy)
grid on
xlabel('t [s]')
ylabel('E [J]')

subplot(2,3,6)
plot(t,mu)
grid on
xlabel('t [s]')
ylabel('\mu [J/T]')

%saveas(fig{i}, output{i}, 'epsc');

end

saveas(fig{1}, 'graphs/gEuler', 'epsc');
saveas(fig{2}, 'graphs/gEulerCromer', 'epsc');
saveas(fig{3}, 'graphs/gRungeKutta2', 'epsc');
saveas(fig{4}, 'graphs/gApp2pos', 'epsc');
saveas(fig{5}, 'graphs/gApp2neg', 'epsc');
