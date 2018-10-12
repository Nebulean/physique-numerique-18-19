#!/usr/bin/env bash

./Exercice2 configuration.in schema=Euler output='Euler.out' nsteps=500
./Exercice2 configuration.in schema=EulerCromer output='EulerCromer.out' nsteps=500
./Exercice2 configuration.in schema=RungeKutta2 output='RungeKutta2.out' nsteps=500
