%% BEFORE-SIMULATIONS
tFin = 172800; % 2 days in seconds.
dtad='true';
atm='true';

nsimul = 20;
%epsilon = logspace(2, -3.5, nsimul);
%epsilon = logspace(2, -4, nsimul);
epsilon = logspace(0, -4, nsimul);

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

amax = zeros(1,nsimul);

tfinamax = zeros(1,nsimul);
for i=1:nsimul
    d = load(output{1,i});
    t = d(:,1);
    ax = d(:,15);
    ay = d(:,16);

    
    a = sqrt(ax.^2 + ay.^2);
    
    [maxAccel, index] = max(a);

    nstepsamax(i) = length(a);

    tfinamax(i) = d(index,1)
  
    if index+2 > length(a)
        fit = polyfit(t(index-4:index), a(index-4:index), 2);
    elseif index-2 < length(a)
        fit = polyfit(t(index:index+4), a(index:index+4), 2);
    else
        fit = polyfit(t(index-2:index+2), a(index-2:index+2), 2);
    end
    A = fit(1); B = fit(2); C = fit(3);

    amax(i) = abs( C - B^2/(4*A) );    
end

amaxref = 2.205352363586426e+02;
erramax = abs(amaxref - amax)

%% GRAPHS
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

xlabel("Number of steps [-]");
ylabel("$|$error on $a_{max}| [m/s^2]$")

set(gca, 'XScale','log');
set(gca, 'YScale','log');

lgd = sprintf("slope: %0.4f", slope);
legend([fitplot], lgd);

box on;
grid on;


hold off;

saveas(fig2, 'graphs/ex2a_conv','epsc');




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
