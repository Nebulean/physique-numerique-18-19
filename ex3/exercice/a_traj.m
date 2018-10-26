set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

d = load("a.out");

t = d(:,1);
theta = d(:,2);
thetadot = d(:,3);

A = 1e-6;
g = 9.81;
L = 0.1;
omega0 = sqrt(g/L)
theta_th = zeros(1,length(t));
for i=1:length(t)
   theta_th = A*cos(omega0*t);
end

f1 = figure
hold on;
p = plot(t, theta, 'Color', [0.8 0.2 0.2], 'LineWidth',1.1);
p_th = plot(t, theta_th, 'Color', [0.2 0.2 0.8], 'LineWidth', 1.1);

set(gca, 'fontsize', 20);

ymin=-1.1e-6; ymax=1.1e-6;
ylim([ymin, ymax]);

legendstr = ["Stormer-Verlet", "Analytical solution"];
legend(legendstr);

xlabel("time $t$ [s]");
ylabel("Angle $\theta$ [rad]") %TODO: C'est bien en radian, non ?

grid on;
hold off;

saveas(f1, 'graphs/a_traj_full','epsc');






f2 = figure
hold on;
p = plot(t, theta, 'Color', [0.8 0.2 0.2], 'LineWidth',1.1);
p_th = plot(t, theta_th, 'Color', [0.2 0.2 0.8], 'LineWidth', 1.1);

set(gca, 'fontsize', 20);

ymin=-1.1e-6; ymax=1.1e-6;
xmin=17; xmax=20;
xlim([xmin, xmax]);
ylim([ymin, ymax]);

legendstr = ["Stormer-Verlet", "Analytical solution"];
legend(legendstr);

xlabel("time $t$ [s]");
ylabel("Angle $\theta$ [rad]") %TODO: C'est bien en radian, non ?

grid on;
hold off;

saveas(f2, 'graphs/a_traj_zoomed','epsc');