# nuclear_physics

This code depends on the [gsl library](https://www.gnu.org/software/gsl/) for integration and interpolation, on the [Ceres library ](http://ceres-solver.org/) to solve non-linear systems of equations and on the [Eigen library](https://eigen.tuxfamily.org/index.php?title=Main_Page) to do some vector/matrix operations Ceres  (it is also a Ceres dependence). `Before running the code, make sure the Makefile points to the correct installation folder of the above libraries.`

This code is being build during my doctorate to calculate various quantities related to nuclear matter in astrophysical conditions, such as:

(1) The equation of state for baryonic matter (nucleons, hyperons and deltas included) in beta-equilibrium. It is possible to include  magnetic fields and/or magnetic moments, but only at T=0. 

(2) The EoS for nucleons w/ fixed proton fraction at any temperatures. Low temperatures (<0.5 MeV) will have convergence problems.

(3) EoS, stability window and phase diagram of density and temperature dependent quark mass model matter and MIT bag model (with and without vector interactions).

(4) Nuclear pasta EoS w/ fixed proton fraction or beta-equilibrium at any temperature.

(5) Electron transport in the pasta: thermal and electric conductivity

If you find any  mistakes or have any doubts, do not hesitate in contacting me: mateusreinke@hotmail.com
