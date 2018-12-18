%% BEFORE-SIMULATIONS
tFin = 172800; % 2 days in seconds.
dtad='true';
atm='true';

nsimul = 20;
epsilon = logspace(2, -6, nsimul);

%% SIMULATIONS
output = cell(1, nsimul);
for i=1:nsimul
     output{1, i} = sprintf("epsilon=%0.10f.out", epsilon(i));
     cmd = sprintf("./Exercice4 configuration.in tFin=%f atm=%s dtad=%s epsilon=%f output=%s", tFin, atm, dtad, epsilon(i), output{1,i});

     disp(cmd);
     system(cmd);
end



%% TRAITEMENT DES DONNEES
nstepsamax = zeros(1,nsimul);
nstepspmax = zeros(1,nsimul);
amax = zeros(1,nsimul);
pmax = zeros(1,nsimul);
for i=1:nsimul
    d = load(output{1,i});
    ax = d(:,15);
    ay = d(:,16);
    power = d(:,17);
    
    % tfin = d(end,1)

    a = sqrt(ax.^2 + ay.^2);
    
    [maxAccel, index] = max(a);
    nstepsamax(i) = length(a(1:index));
   
    [maxPower, index] = max(abs(power));
    nstepspmax(i) = length(power(1:index));
    
    amax(i) = maxAccel;
    pmax(i) = maxPower;
end

[tmp, index] = max(nstepsamax);
amaxref = amax(index);
erramax = abs(amaxref - amax)

[tmp, index] = max(nstepspmax(i));
pmaxref = pmax(index);
errpmax = abs(pmaxref - pmax);

%% GRAPHS
fig1=figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(nstepsamax, amax, 'x');

% set(gca, 'XScale','log');
% set(gca, 'YScale','log');

box on;
grid on;


hold off;





fig2=figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(nstepsamax(1:end-1), erramax(1:end-1), 'x');

[fit, slope] = poly_approx(log10(nstepsamax(1:end-1)), log10(erramax(1:end-1)), 1, 2);
fit = 10.^fit;
fitplot = plot(fit(1,:), fit(2,:), '-');


set(gca, 'XScale','log');
set(gca, 'YScale','log');

lgd = sprintf("slope: %0.4f", slope);
legend([fitplot], lgd);

box on;
grid on;


hold off;



fig3=figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(nstepspmax, pmax, 'x');

% set(gca, 'XScale','log');
% set(gca, 'YScale','log');

box on;
grid on;


hold off;




fig4=figure;
hold on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);
set(gca, 'fontsize', 22);

plot(nstepspmax(2:end), errpmax(2:end), 'x');

[fit, slope] = poly_approx(log10(nstepspmax(2:end)), log10(errpmax(2:end)), 1, 2);
fit = 10.^fit;
fitplot = plot(fit(1,:), fit(2,:), '-');


set(gca, 'XScale','log');
set(gca, 'YScale','log');

lgd = sprintf("slope: %0.4f", slope);
legend([fitplot], lgd);

box on;
grid on;


hold off;





function [polynome, slope] = poly_approx(x, y, ordre, steps)
    pf = polyfit(x, y, ordre);
    slope = pf(1);
    T = linspace(min(x), max(x), steps);

    n = ordre + 1;

    polynome = zeros(2,length(T));
    for i=1:n
       polynome(2,:) = polynome(2,:) + pf(i)*T.^(n-i);
    end

    polynome(1,:) = T;
end
