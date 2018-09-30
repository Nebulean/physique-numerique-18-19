% On charge le fichier
d = load("C_all.out");

% On trie les colonnes
%tfin = d(:,1);
zfin = d(:,2);
%vfin = d(:,3);
vinit = d(:,4);

% On calcul la médiane
zMedian = median(zfin);

% On créé un nouveau vecteur, en ne copiant que les valeurs sous la
% médiane.

zCell = [];
vCell = [];

for i=1:length(zfin)
    if zfin(i) < zMedian
       zCell{end+1,1} = zfin(i);
       vCell{end+1,1} = vinit(i);
    end
end

% Maintenant on rempli deux array avec ces valeurs.
z = double.empty(length(zCell),0);
v = double.empty(length(vCell),0);

for i=1:length(zCell)
    z(i) = zCell{i};
    v(i) = vCell{i};
end

f = figure;
hold on;
% On dessine le graph
plot(v, z, '-');
xlabel("vitesse initiale");
ylabel("position finale");

% On ajoute une ligne horizontale où se trouve la lune
refline(0, 384400000)
hold off;