tFin = 172800; % 2 days in seconds.
dtad="true";
atm="false";
dt=60; % dt initial

epsilon = 1e-5;

cmd = sprintf("./Exercice4 configuration.in output=dtevolution.out tFin=%f dtad=%s dt=%f epsilon=%f atm=%s", tFin, dtad, dt, epsilon, atm);
system(cmd);

d = load("dtevolution.out");

t = d(:,1);
dt = d(:,2);

fig1=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(t, dt, '.');
grid on;
box on;

ylabel("Time step $\Delta t$ [s]");
xlabel("Time $t$ [s]");

hold off;

saveas(fig1, 'graphs/ex1c_conveps_dt','epsc');

nsteps = length(t)