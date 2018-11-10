d = load("d_thm.out");

% t theta thetadot emec pnc emecdot
t = d(:,1);
theta = d(:,2);
thetadot = d(:,3);
emec = d(:,4);
pnc = d(:,5);
emecdot = d(:,6);

f1 = figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);


emecdot2 = diff(emec)./diff(t);
newpnc = pnc;
newpnc(end) = [];
newt = t;
newt(end) = [];

plot(newt, abs(emecdot2-newpnc));

%ymin=-2e-16;
%ymax=2e-16;
%ylim([ymin, ymax]);

xlabel("t [s]");
ylabel("$\left|\frac{d(E_m)}{dt} - P_{nc}\right|$ [W]")

grid on;
hold off;

f2 = figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(t,theta);

xlabel("t [s]");
ylabel("$\theta$ [rad]")
grid on;

hold off;

%% graph emecdot
f3 = figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

%emecdot2 = diff(emec)./diff(t);
%t(end) = [];
plot(t,emecdot);

xlabel("t [s]");
ylabel("derivative of mecanical energy [J]")
grid on;

hold off;

saveas(f1, 'graphs/d_thm','epsc');
saveas(f2, 'graphs/d_theta','epsc');
