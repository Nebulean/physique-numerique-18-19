d = load("a.out");

t = d(:,1);
dt = d(:,2);

x1 = d(:,3);
y1 = d(:,4);

x2 = d(:,7);
y2 = d(:,8);

x3 = d(:,11);
y3 = d(:,12);

%% On plot les endroits initiaux.
fig=figure
hold on;
plot(x1(1), y1(1), 'o', 'Color','blue');
plot(x2(1), y2(1), 'o', 'Color','red');
plot(x3(1), y3(1), 'o', 'Color','green');

% puis on plot les positions.
plot(x1(2:end), y1(2:end), '.', 'Color', 'blue');
plot(x2(2:end), y2(2:end), '.', 'Color', 'red');
plot(x3(2:end), y3(2:end), '.', 'Color', 'green');

hold off;