#!/usr/bin/env bash

echo "Simulation de l'exercice a."
#a
./Exercice3 configuration.in kappa=0. d=0. theta0=1e-6 thetadot0=0. tfin=20 output=a_traj.out
./Exercice3 configuration.in kappa=0. d=0. theta0=1e-6 thetadot0=0. tfin=20 output=a_ener.out

echo "Simulation de l'exercice B."

#b
./Exercice3 configuration.in kappa=0. d=0. thetadot0=0 tFin=20 output=b.out

#c
./Exercice3 configuration.in theta0=0. thetadot0=1e-2 Omega=9.9045 d=0.03 kappa=0. tFin=250 output=c_thm.out

#d
./Exercice3 configuration.in theta0=0. thetadot0=1e-2 Omega=19.809 d=0.005 kappa=0.05 tFin=100 output=d_thm.out

#merde
./Exercice3 configuration.in theta0=3.14 thetadot0=1e-2 Omega=9.9045 d=0.04 kappa=0. tFin=100 output=merde.out

#g
./Exercice3 configuration.in theta0=1.0471975512 thetadot0=1e-2 Omega=9.9045 d=0.05 kappa=0.03 tFin=100 output=g_fourier.out
# pi/3  = 1.0471975512
# pi    = 3.1415926536
