%% On génère les données
dtad="false"; % Pour ex1b
%dtad="true"; % Pour ex1c
cmd = sprintf("./Exercice4 configuration.in tFin=172800 atm=false dt=8 dtad=%s output=ex1b_traj.out", dtad);
system(cmd)
%% On load les données
d = load("ex1b_traj.out");

t = d(:,1);
dt = d(:,2);

x1 = d(:,3);
y1 = d(:,4);

x2 = d(:,7);
y2 = d(:,8);

x3 = d(:,11);
y3 = d(:,12);

vx3 = d(:,13);
vy3 = d(:,14);

%% On plot les endroits initiaux.
fig=figure
hold on;
pbaspect([1 1 1]);
daspect([1 1 1])
plot(x3(1), y3(1), 'x', 'Color','green');
centerOfEarth = [x1(1), y1(1)];
plotCircle(centerOfEarth, 6371000, 500, 'blue');

% puis on plot les positions.
plot(x1(2:end), y1(2:end), '-', 'Color', 'blue', 'LineWidth', 1.2);
plot(x3(2:end), y3(2:end), '-', 'Color', 'green', 'LineWidth', 1.2);

grid on;


hold off;

nsteps = length(t)

function circle = plotCircle(center, radius, nb, color)
    circle = zeros(nb, 2);
    t = linspace(0, 2*pi, nb)

    for i=1:nb
        circle(i,1) = center(1) + radius*cos(t(i));   %*(1-t.^2)/(1+t.^2);
        circle(i,2) = center(2) + radius*sin(t(i));   %*2.*t/(1+t.^2);
    end
       
    plot(circle(:,1), circle(:,2), '-', 'Color', color, 'LineWidth', 1.2);
end


% function polynome = poly_approx(x, y, ordre, steps)
%     pf = polyfit(x, y, ordre);
%     T = linspace(min(x), max(x), steps);
%     
%     n = ordre + 1
%     
%     polynome = zeros(2,length(T))
%     for i=1:n
%        polynome(2,:) = polynome(2,:) + pf(i)*T.^(n-i);
%     end
%     polynome(1,:) = T;
% end