d = load("2corps.out");

t = d(:,1);
dt = d(:,2);

e1 = d(:,18);
e2 = d(:,19);

p1 = d(:,20);
p2 = d(:,21);

dist = d(:,22);

f1=figure
plot(t,e1);
grid on;

f2=figure
plot(t,e2);
grid on;

f3=figure
plot(t,p1);
grid on;

f4=figure
plot(t,p2);
grid on;

f5=figure
plot(t,dist);
grid on;