# BSPCI2e (B-splines Configuration Interaction 2-electron) -            The 2-electron system simulator
#
Overview: This code is capable of simulating the system of two interacting electron that are trapped in some form of quantum well (Nuclear, Quantum-Dot, etc) and interacting with a high intensity and ultrashort laser pulse in the XUV wavelength regime. 

This is achieved by first solving the single electron field-free (no laser) Schrodinger equation by converting it to a matrix equation using a basis of B-spline functions and applying a diagonalization routine to find the energies, coefficients and radial wave-functions. These functions become the function basis for the two electron field free schrodinger equation which is again diagonalized on a basis of configurations (products of one-electron basis functions) to retreive the two-electron energies, configuration- interaction coefficients and two electron wavefunctions. Finally these two electron functions become the basis for the time dependent Shrodinger equation which includes the full system (two interacting electrons + laser). This is solved by projection conversion into a massive system of ordinary differential equations to retrieve time dependent coefficients. All of the coefficients and energies retrieved up to this point can be used to generate obervable properties of the system such as radial distributions, angular distributions, and energy distributions, to name a few.


Required libraries to link:

LAPACK - Linear Algebra PACKage, 

BLAS - Basic Linear Algebra Subroutines,

NAG - Numerical Algorithms Group,

NetCDF - Network Common Data Form

Some parameters this code recieves are physical parameters of the laser i.e photon energy, laser intensity, and laser pulse duration. Once the libraries are linked and the configuration files (basis functions) are set up, a script can be run to solve the electronic structure of the field free system, and another script to run the laser pulse propagation step. Various other programmes evaluate the overbales such as the distributions mentioned above and expectation values, etc  

Various inputs into this code are non-physical parameters such as the 'Size of the box' that the solution is solved in. In real life. the system lives in the universe which extends effectively to infinity, obviously a computer cannot simulate this so we choose a box size large enough that the effect of the edge of the box is negligible on the system.

Other parameters that affect the solution are the number of basis functions we choose. We might choose 60 basis-spline functions to generate 60 single electron solutions on a single CPU. The more we choose, the more accurate our solutions. However the final steps of this method (configuration interaction method) are highly computationally intensive. Even this tiny 2-electron system is enough to cause a non-paralellized code of a very fast desktop computer to take several hours to evalute the atomic structure. So we set limits on how many functions we use and the size of the box.
