%% Préparation de(s) simulation(s)
schema = ["A", "B", "C"];
omega = 2*pi/900;
n_stride = 10;
Npoints = 5000;
tfin = 10000;
pulse = "true";

for i=1:3
    param(i) = sprintf("schema=%s n_stride=%i Npoints=%i tfin=%i omega=%f pulse=%s", schema(i), n_stride, Npoints, tfin, omega, pulse);
    output(i) = sprintf("tsunami_a_%s", schema(i));
end


set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);


%% Simulations
for i=1:3
    cmd = sprintf("./Exercice7 configuration_tsunami.in output=%s %s", output(i), param(i));
    disp(cmd);
    system(cmd);
end

x1wave = {};
t1wave = {};
u_calc = {};
for j=1:3
    %% Chargement des données
    dataf = load(sprintf("%s_f.out", output(j)));
    datau = load(sprintf("%s_u.out", output(j)));
    dataE = load(sprintf("%s_E.out", output(j)));

    %% Traitement des données
    f = dataf(:,2:end);
    t = dataf(:,1);

    x = datau(:,1);

    % On veut suivre une unique vague, du coup, je dois trouver le maximum de
    % la première vague. On selectionne le dernier pic pour chaque instant.
    x1wave{j} = [];
    t1wave{j} = [];

    % Méthode 1 - maximum locaux (mode pulse et continu)
    xbefore = 0;
    for i=1:length(f(:,1))
        indexes = [];
        % On trouve les peaks à un instant
        [pks, indexes] = findpeaks(f(i,:));


        if ~isempty(indexes)
            disp("Is not empty");
            % On check si y'a eu retour en arrière
            if (xbefore > x(indexes(end)))
               disp("Stopping for");
               fprintf("Current x = %f", x(indexes(end)));
               break;
            else
                disp("Continuing for");
                xbefore = x(indexes(end));
            end

            % On note la position à l'instant t.
            x1wave{j}(i) = x(indexes(end)); 
            t1wave{j}(i) = t(i);
        else
            disp("Is empty");
        end
    end

    % Méthode 2 - maximum (mode pulse) -> Marche pas.
%     for i=1:length(f(:,1))
%         [M, index] = max(f(i,:));
%         x1wave{j}(i) = x(index);
%         t1wave{j}(i) = t(i);
%     end
    

    u_calc{j} = diff(x1wave{j})./diff(t1wave{j});
end


% On créé la valeur théorique
g=9.81;
uth = sqrt(g*h(x));

%% On fait des graphes
fig_u=figure
hold on;

set(gca, 'fontsize', 25);
set(gca, 'LineWidth',1.5);

plot(x(1:end-1), uth, '-', 'LineWidth', 1.5); % Vitesse théorique
for j=1:3
    plot(x1wave{j}(1:end-1), u_calc{j}, '-.', 'LineWidth', 1.5); % Vitesse expérimentale
end

xlabel("$x~[m]$");
ylabel("$u~[m/s]$");

legend(["th", "A", "B", "C"]);

%plot(x, f(indexes(end),:))
% plot(x(indexes), f(end, indexes),'o')

grid on;
box on;

hold off;







function res = h(x)
    xa = 200000;
    xb = 370000;
    xc = 430000;
    xd = 600000;
    L = 800000;
    
    hrecif = 20;
    hocean = 8000;
    
    res = [];
    
    for i=1:length(x)
        if (0 <= x(i)) & (x(i) < xa)
            res(i) = hocean;
        elseif (xa <= x(i)) & (x(i) < xb)
            res(i) = hocean + (hrecif-hocean) * sin(pi*(x(i)-xa)/(2*(xb-xa)))^2;
        elseif (xb <= x(i)) & (x(i) < xc)
            res(i) = hrecif;
        elseif (xc <= x(i)) & (x(i) < xd)
            res(i) = hrecif - (hrecif - hocean) * sin(pi*(xc-x(i))/(2*(xc-xd)))^2;
        elseif (xd <= x(i)) & (x(i) <= L)
            res(i) = hocean;
        end
    end
end