d = load("a.out");

t = d(:,1);
dt = d(:,2);

x1 = d(:,3);
y1 = d(:,4);

% x2 = d(:,7);
% y2 = d(:,8);

x3 = d(:,11);
y3 = d(:,12);

vx3 = d(:,13);
vy3 = d(:,14);
%% On plot les endroits initiaux.
fig=figure
hold on;
plot(x3(1), y3(1), 'x', 'Color','green');
centerOfEarth = [x1(1), y1(1)];
plotCircle(centerOfEarth, 6371000, 50, 'blue');

centerOfMoon = [x2(1), y2(1)];
plotCircle(centerOfMoon, 1737500, 50, 'red');


% puis on plot les positions.
plot(x1(2:end), y1(2:end), '-', 'Color', 'blue', 'LineWidth', 1.2);
% plot(x2(2:end), y2(2:end), '-', 'Color', 'red');
plot(x3(2:end), y3(2:end), '-', 'Color', 'green', 'LineWidth', 1.2);



hold off;

function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb)
    
    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
   
%     circle(:,1) = center(1) + radius.*cos(t)   %*(1-t.^2)/(1+t.^2);
%     circle(:,2) = center(2) + radius.*sin(t)   %*2.*t/(1+t.^2);
    
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end