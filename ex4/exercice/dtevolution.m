d = load("a.out");

t = d(:,1);
dt = d(:,2);

fig1=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(t, dt, 'x');
grid on;

ylabel("$\Delta t$ [s]");
xlabel("Time $t$ [s]");