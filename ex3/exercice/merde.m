m = load("merde.out");

t = m(:,1);
theta = m(:,2);
thetadot = m(:,3);
emec = m(:,4);
pnc = m(:,5);
emecdot = m(:,6);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18); 

plot(t,theta)