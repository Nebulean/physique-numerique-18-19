%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = 'output';
data = load([fichier,'_u.out']);
x = data(:,1);
u = data(:,2);
data = load([fichier,'_E.out']);
t = data(:,1);
E = data(:,2);
data = load([fichier,'_f.out']);
f = data(:,2:end);

%% Figures %%
%%%%%%%%%%%%%
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

% figure('Name',['Analyse de ' fichier])
% subplot(2,2,1);
figu=figure;
plot(x,u)
grid
xlabel('$x$ [m]')
ylabel('$u$ [m/s]')

% subplot(2,2,2);
figE=figure;
plot(t,E)
grid
xlabel('$t$ [s]')
ylabel('$E$ [m$^3$]')

% xt = subplot(2,2,4);
figxt=figure;
pcolor(x,t,f)
shading interp
colormap jet
c = colorbar;
xlabel('$x$ [m]')
ylabel('$t$ [s]')
ylabel(c,'f(x,t) [m]')

% fct = subplot(2,2,3);
figf=figure;
h = plot(x,f(1,:));
grid
xlabel('$x$ [m]')
ylabel('$f(x,t)$ [m]')
ht = title('$t=0$ s');
ylim([min(f(:)),max(f(:))])
for i=2:length(t)
    pause(.05)
    if ~ishandle(h)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(h,'YData',f(i,:))
    set(ht,'String',sprintf('$t=%0.2f$ s',t(i)))

end

saveas(figE,'graphs/ex1Efixe','epsc');
saveas(figxt,'graphs/ex1xtfixe','epsc');
saveas(figf,'graphs/ex1ffixe','epsc');

% saveas(figE,'graphs/ex1Elibre','epsc');
% saveas(figxt,'graphs/ex1xtlibre','epsc');
% saveas(figf,'graphs/ex1flibre','epsc');