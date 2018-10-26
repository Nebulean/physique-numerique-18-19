#!/usr/bin/env bash

echo "Simulation de l'exercice a."
#a
./Exercice3 configuration.in kappa=0. d=0. theta0=1e-6 thetadot0=0. tfin=20 output=a.out

echo "Simulation de l'exercice B."

#b
./Exercice3 configuration.in kappa=0. d=0. thetadot0=0 tFin=20 output=b.out
