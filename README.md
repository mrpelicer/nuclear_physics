# nuclear_physics

This code depends on the [gsl](https://www.gnu.org/software/gsl/) and [Cuba](https://feynarts.de/cuba/) libraries for integration and interpolation, on the [Ceres](http://ceres-solver.org/) library to solve non-linear systems of equations and on the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library to do some vector/matrix operations in the pasta code. 

The code is tested was in Ubuntu 14.04, 16.04 and 22.04. On them you can install the dependencies (except cuba) by typing:

`sudo apt install libgsl-dev libeigen3-dev libceres-dev`

For the Cuba library, follow the instructions in their page.

`Make sure the Makefile points to the correct installation folder of the above libraries if there is a problem with the compilation`


This code is being build during my doctorate to calculate various quantities related to nuclear matter in astrophysical conditions, such as:

(1) The equation of state (EoS) for baryonic matter (with nucleons, hyperons and/or deltas included) in beta-equilibrium. It is possible to include  magnetic fields and/or magnetic moments, but only at T=0. 

(2) The EoS for nucleons w/ fixed proton fraction at any temperatures. Very low temperatures (<0.1 MeV) will have convergence problems: you will have to tackle with the gsl parameters of integration.

(3) EoS, stability window and phase diagram of density and temperature dependent quark mass model matter and vector MIT bag model.

(4) Nuclear pasta properties w/ fixed proton fraction or beta-equilibrium at any temperature.

(5) Electron transport in the pasta: thermal and electric conductivities.

If you find any  mistakes or have any doubts, do not hesitate in contacting me: mateusreinke@hotmail.com

In case you find the code useful, consider citing my publications related to your work: see [INSPIRE](https://inspirehep.net/authors/1905850)

If you find any  mistakes or have any doubts, do not hesitate in contacting me: mateusreinke@hotmail.com
