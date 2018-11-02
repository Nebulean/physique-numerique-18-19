c = load("c_thm.out");

t = c(:,1);
theta = c(:,2);
thetadot = c(:,3);
emec = c(:,4);
pnc = c(:,5);
emecdot = c(:,6);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

f1 = figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);
plot(t, emecdot-pnc);
xlabel("t [s]");
ylabel("$\frac{dP}{dt} - P_{nc}$ [W]")
grid on;
hold off;

f2 = figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
plot(t,theta);
xlabel("t [s]");
ylabel("$\theta [rad]$")
set(gca, 'fontsize', 22);
grid on;
hold off;

saveas(f1, 'graphs/c_thm', 'epsc');
saveas(f2, 'graphs/c_theta', 'epsc');