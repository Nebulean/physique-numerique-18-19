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
emecdot2 = diff(emec)./diff(t);
newpnc = pnc;
newpnc(end) = [];
newt = t;
newt(end) = [];
plot(newt, abs(emecdot2-newpnc));
set(gca, 'fontsize', 22);
xlabel("t [s]");
ylabel("$\frac{dE_m}{dt} - P_{nc}$ [W]")
grid on;
hold off;

f2 = figure;
hold on;
plot(t,theta);
set(gca, 'fontsize', 22);
xlabel("t [s]");
ylabel("$\theta [rad]$")
grid on;
hold off;

%% graph emecdot
f3 = figure;
hold on;
plot(newt,emecdot2);
set(gca, 'fontsize', 22);
xlabel("t [s]");
ylabel("derivative of mechanical energy [J]")
grid on;

hold off;

%% graph Pnc
f4 = figure;
hold on
plot(t,pnc);
set(gca, 'fontsize', 22);
xlabel("t [s]");
ylabel("$P_{nc}$")
grid on
hold off

%% graph emec
f4=figure;
plot(t,emec)
grid on
set(gca, 'fontsize', 22)
xlabel("t [s]")
ylabel("$E_{mec}$")

%%
saveas(f1, 'graphs/c_thm', 'epsc');
saveas(f2, 'graphs/c_theta', 'epsc');