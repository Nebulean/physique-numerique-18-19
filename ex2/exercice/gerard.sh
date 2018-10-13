#!/usr/bin/env bash

./Exercice2 configuration.in schema=Euler output='Euler.out' nsteps=500
./Exercice2 configuration.in schema=EulerCromer output='EulerCromer.out' nsteps=500
./Exercice2 configuration.in schema=RungeKutta2 output='RungeKutta2.out' nsteps=500

# Pour l'application 1
# ii)
./Exercice2 configuration.in schema=RungeKutta2 output='app1_ii.out' nsteps=10000

# Pour l'application 2
# i)
./Exercice2 configuration.in Kappa=100. schema=RungeKutta2 output='App2pos.out' nsteps=10000
./Exercice2 configuration.in q=-1.6022e-19 Kappa=100. x0=0.00139192 schema=RungeKutta2 output='App2neg.out' nsteps=10000
