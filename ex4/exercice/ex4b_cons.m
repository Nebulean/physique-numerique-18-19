d = load("2corps.out");

t = d(:,1);
dt = d(:,2);

e1 = d(:,18);
e2 = d(:,19);

p1 = d(:,20);
p2 = d(:,21);

dist = d(:,22);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

f1=figure
plot(t,e1);
xlabel('Time [s]');
ylabel('Mechanical Energy [J]');
set(gca, 'fontsize', 22);
box on;
grid on;

f2=figure
plot(t,e2);
% legend('Earth','Moon');
box on;
xlabel('Time [s]');
ylabel('Mechanical Energy [J]');
set(gca, 'fontsize', 22);
grid on;

f3=figure
plot(t,p1, t,p2);
legend('Earth','Moon');
box on;
xlabel('Time [s]');
ylabel('Momentum [kg m/s]');
set(gca, 'fontsize', 22);
grid on;

f4=figure
plot(t,dist);
box on;
xlabel('Time [s]');
ylabel('Distance [m]');
set(gca, 'fontsize', 22);
grid on;

saveas(f1,'graphs/ex4b_emece.eps','epsc');
saveas(f2,'graphs/ex4b_emecl.eps','epsc');
saveas(f3,'graphs/ex4b_p.eps','epsc');
saveas(f4,'graphs/ex4b_d.eps','epsc');