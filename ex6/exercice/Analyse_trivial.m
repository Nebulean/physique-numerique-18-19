function Analyse_trivial

filename='output';

data = load([filename '_Er_Dr.out']);
rmid = data(:,1);
Er = data(:,2);
Dr = data(:,3);
data = load([filename '_phi.out']);
r = data(:,1);
phi = data(:,2);
data = load([filename '_rholib_divEr_divDr.out']);
rmidmid = data(:,1);
rholib = data(:,2);
divEr = data(:,3);
divDr = data(:,4);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

%Cas trivial
R=.12;
phith=(R^2-r.^2)/4;

f1=figure;
plot(r,phi,r,phith, 'LineWidth',1.5)
xlabel('$r$ [m]')
ylabel('$\phi$ [V]')
grid on, box on
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);
legend('numerical result','analytical result');

f2=figure;
hold on
plot(rmid,Er, 'LineWidth',1.5)
xlabel('$r$ [m]')
ylabel('$E_r$ [V/m]')
grid on, box on
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

f3=figure;
hold on
plot(rmid,Dr, 'LineWidth',1.5)
xlabel('$r$ [m]')
ylabel('$D_r/\epsilon_0$ [V/m]')
grid on, box on
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

f4=figure;
hold on
plot(rmidmid,rholib,'DisplayName','$\rho_{lib}/\epsilon_0$', 'LineWidth',1.5)
plot(rmidmid,divDr,'--','DisplayName','$div(D_r)/\epsilon_0$', 'LineWidth',1.5)
plot(rmidmid,divEr-divDr,'DisplayName','$\rho_{pol}/\epsilon_0$', 'LineWidth',1.5)
xlabel('$r$')
ylabel('$\rho/\epsilon_0$ [V/m$^2$]')
legend('show')
grid on, box on
set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

saveas(f1, "graphs/c_phi","epsc");
saveas(f2, "graphs/c_Er","epsc");
saveas(f3, "graphs/c_Dr","epsc");
% saveas(f4, "graphs/c_rho_div","epsc");
