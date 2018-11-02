d = load("d_thm.out");

t = d(:,1);
theta = d(:,2);
thetadot = d(:,3);
emec = d(:,4);
pnc = d(:,5);
emecdot = d(:,6);

f = figure;
hold on;
plot(t, emecdot-pnc);

xlabel("t [s]");
ylabel("$\frac{dP}{dt} - P_{nc}$ [W]")

grid on;
hold off;