clear all;
format long;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);

%% Initialisation
valb=0.02;
valR=0.12;

nsimul=20;
N=round(logspace(2, 5, nsimul));
valN1=N;
valN2=2.*N;

param = {};
cmd = {};
output = {};

%% Simulation
for i=1:nsimul
    output{i} = sprintf("outputdii_%i", i);
    param{i} = sprintf("output=%s trivial=false b=%f R=%f N1=%i N2=%i",output{i} ,valb,valR,valN1(i),valN2(i));
    cmd{i} = sprintf("./Exercice6 configuration.in %s",param{i});
    disp(cmd{i});
    %system(cmd{i});
end


%% On importe et traite les données
for i=1:nsimul
    clear data M index rmidmid divEr divDr;
    loading = sprintf("%s_rholib_divEr_divDr.out", output{i})
    data = load(loading);
    rmidmid = data(:,1);
    divEr = data(:,3);
    divDr = data(:,4);
    
    [M, index] = max(abs(divEr - divDr));

    chargePol(i) = 8.85418782e-12*((divEr(index) - divDr(index)) * (rmidmid(index + 1) - rmidmid(index - 1)))/2; % approximation par un triangle 
end

%% On plot le résultat.
figPolCh=figure;
hold on;

plot(1./N, chargePol, '+', 'MarkerSize', 10);

fit = poly_approx(1./N, chargePol, 1, 2);
plot(fit(:,1), fit(:,2), '-', 'LineWidth', 1.2);

box on;
grid on;

xlabel("1/N");
ylabel("Polarization charge in $r=b$ [$C$]");

legend(["Data", "Linear fit"], 'Location', 'southeast');

set(gca, 'LineWidth',1.5);
set(gca, 'fontsize',25);

hold off;

saveas(figPolCh, "graphs/exdii-conv-polch", "epsc");

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
    
    polynome=polynome.'; % transposition pour améliorer l'utilisation.
end