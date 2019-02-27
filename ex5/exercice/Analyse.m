%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = 'output';
data = load([output '_T.out']);
N = sqrt(length(data));
Y = data(1:N,2);
X = data(1:N:N*N,1);
T = reshape(data(:,3),N,N)'; % 1D -> 2D
h = (max(X)-min(X))/(N-1);

data = load([output '_P.out']);
t = data(:,1);
Pc = data(:,2);
Pf = data(:,3);
Ptot = data(:,4);
kappa=1.2;

data = load(['output_F.out']);
% Yc = data(1:N,2);
% Xc = data(1:N:N*N,1);
jyc = reshape(data(:,3),N-1,N-1)';
jxc = reshape(data(:,4),N-1,N-1)';
% jxc2 = reshape(data(:,5),N-1,N-1)';
% jyc2 = reshape(data(:,6),N-1,N-1)';


%% Analyse %%
%%%%%%%%%%%%%
% TODO : calcul du flux de chaleur au centre des cellules du maillage
Xmid = X(1:N-1)+h/2;
Ymid = Y(1:N-1)+h/2;
% jxc = zeros(N-1,N-1);
% jyc = zeros(N-1,N-1);
% jx = zeros(N-1,N-1);
% jy = zeros(N-1,N-1);

% Simple différence finie centrée, mais on veut la valeur au centre de la
% cellule.
% for i=1:length(jx)
%    for j=1:length(jx)
%       jxc(i,j) = -kappa/h * (T(i+1,j) - T(i,j));
%       jyc(i,j) = -kappa/h * (T(i,j+1) - T(i,j));
%    end
% end

% for i=1:length(jxc)
%    for i=1:length(jyc)
%       jxc(i,j) = ;
%       jyc(i,j) = ;
%    end
% end

% for i=2:N-1
%     for j=2:N-2
%         jx(i-1,j) = -kappa/h * (T(i+1,j) - T(i,j));
%         jy(i-1,j) = -kappa/h * (T(i,j+1) - T(i,j));
%     end
% end

% Du coup, on centre le résultat.
% for i=1:N-1
%     for j=2:N-2
%         jyc(i,j) = (jxc(i,j+1)+jx(i,j))/2
%         jxc(i,j) = (jx(i+1,j)+jx(i,j))/2;
%     end
% end
jnorm = sqrt(jxc.^2+jyc.^2);

%% Figures %%
%%%%%%%%%%%%%
% Temperature :
f1=figure
hold on
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'defaultLegendInterpreter', 'latex');
% set(groot, 'defaultTextInterpreter', 'latex');
% set(groot, 'defaultAxesFontSize', 18);
% set(gca, 'fontsize', 30);
% set(gca, 'LineWidth',1.5);

contourf(X,Y,T',15,'LineStyle','None'), hold on
stride = 2; % (affiche 1 point sur [stride] dans chaque dimension)
quiver(Xmid(1:stride:end,1:stride:end),Ymid(1:stride:end,1:stride:end),jxc(1:stride:end,1:stride:end)',jyc(1:stride:end,1:stride:end)','k')
xlabel('x [m]')
ylabel('y [m]')
title('T(x,y) [°C]')
colorbar
axis equal
box on
hold off

% Flux de chaleur :
f2=figure
hold on


contourf(Xmid,Ymid,jnorm',30,'LineStyle','None')
xlabel('x [m]')
ylabel('y [m]')
title('|j|(x,y) [W/m]')
colorbar
axis equal
box on
hold off


% Puissance :
f3=figure
hold on
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'defaultLegendInterpreter', 'latex');
% set(groot, 'defaultTextInterpreter', 'latex');
% set(groot, 'defaultAxesFontSize', 18);
% set(gca, 'fontsize', 30);
% set(gca, 'LineWidth',1.5);

plot(t, Pc, t, Pf, t, Pc + Pf, t, Ptot)
xlabel('t [s]')
ylabel('P [W]')
legend('P_c', 'P_f', 'P_c+P_f', 'P_{tot}')
grid on
box on
hold off

saveas(f1, "graphs/c_temp","epsc");
saveas(f2, "graphs/c_heat_flux","epsc");
saveas(f3, "graphs/c_power","epsc");
