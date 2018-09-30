rm C_all.out
touch C_all.out

# for i in {11000..11200}
for (( i=11000; c<=11200; i+=1))
do
  # On commence par la simulation
  ./Exercice1 configurationA.in v0=$i nsteps=64000 rho0=0 output="tmp.out" tfin=350400

  # On obtient la dernière valeur ajoutée
  LAST=$(tail -1 tmp.out)

  # On ajoute la vitesse actuelle
  LAST="$LAST $i"

  # On vérifie si la string contient "nan"
  if [[ ${LAST} != *nan* ]];
  then
    # echo $LAST
    echo $LAST >> C_all.out
  fi
  echo $LAST
  # echo "$i"
  # echo $LAST
done

# 97h20min = 350400s
