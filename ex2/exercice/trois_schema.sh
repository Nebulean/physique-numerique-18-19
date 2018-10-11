#!/usr/bin/env bash

./Exercice2 configuration.in schema=Euler output='Euler.out' steps=10000
./Exercice2 configuration.in schema=EulerCromer output='EulerCromer.out' steps=10000
./Exercice2 configuration.in schema=RungeKutta2 output='RungeKutta2.out' steps=10000
