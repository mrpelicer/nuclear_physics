# nuclear_physics

This code depends on the [gsl library](https://www.gnu.org/software/gsl/) for integration and interpolation, on the [Ceres library ](http://ceres-solver.org/) to solve systems of equations, and on the [Eigen library]() to do some matrix operations Ceres and as a Ceres dependence. `Make sure that the Makefile points to the correct installation folders in your computer.`

This code is being build during my doctorate to calculate various quantities, such as:

(1) The equation of state for baryons (nucleons, hyperons, deltas)in beta-equilibrium. It is possible to include  magnetic fields and/or magnetic moments, but only at T=0 for now. 

(2) The EoS for nucleons w/ fixed proton fraction at any temperature temperature.

(3) EoS and radius of the nuclear pasta, w/ fixed proton fraction or beta-equilibrium at any temperature.

(4) EoS, stability window and phase diagram of density and temperature dependent quark mass model matter.

If you find any  mistakes, do not hesitate in contacting me: mateusreinke@hotmail.com