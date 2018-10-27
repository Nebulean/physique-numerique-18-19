set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

d = load("a_ener.out");

t = d(:,1);
ene = d(:,4);

f = figure;
hold on;
plot(t, ene, '-', 'LineWidth', 1.2);

xlabel('Time t [s]');
ylabel('Energy $E$ [J]');

set(gca, 'fontsize', 20);

ymin = -.098099999999955;
ymax = -.098099999999950;
ylim([ymin, ymax])

grid on;

hold off;

saveas(f, 'graphs/a_ener','epsc');