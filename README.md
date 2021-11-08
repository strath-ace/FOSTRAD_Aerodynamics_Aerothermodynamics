# FOSTRAD

(The Free Open Source Tool for Re-entry of Asteroids and Debris)

`FOSTRAD` is a simulation suite that allows for the estimation of aerodynamics and aerothermodynamics of an object experiencing flow-fields during atmospheric re-entry scenarios. It has better prediction capability than object-oriented re-entry codes as it employs local panel formulation to compute aerodynamic and aerothermodynamic coefficients. The unique capabilities of the code are described below.

#### 1. Grid Refinement
An internal mesh refiner is used to ensure that there are enough facets to provide accurate aerodynamic and aerothermodynamic coefficients.

#### 2. Visible Facet Detection
The visibility facet detection algorithm employs a fast and robust formulation based on occlusion culling and back-face culling. This approach gives the capability to detect the shadowed facets based on object orientation and self-shadowing of facets.

#### 3. Local radius Computation
A local radius estimation algorithm is implemented to provide better and consistent heat-flux computations.

#### 4. Shape-based Correction Factor
Area-ratio based drag and pitching moment correction factors are applied for cylindrical and parallelepiped geometries based on DSMC simulations.

#### 5. Aerothermodynamic Models,
Stagnation point heat transfer in the continuum flow regime can be computed using four different semi-empirical models.
1. SCARAB model (`sc`) : Provides higher heat flux for TPS based spacecraft.
2. Detra-Kemp-Riddel model (`krd`): Uses a Reynolds number formulation with effective nose radius of the object and an assumption of super catalytic wall boundary condition.
3. Fay-Riddel model (`fr`): Includes dissociation effects with fully catalytic boundary conditions.
4. Van-Driest model (`vd`): Simplified model of Fay-Riddel formulation.

All the above models use a local inclination angle to compute overall heat transfer.

#### 6. Transitional Bridging Function
A local radius based bridging model as function of Knudsen number and effective nose radius is implemented. This model has been calibrated to DSMC results to improve accuracy.


Requirements
------

`FOSTRAD` is currently written in `MATLAB` and hence needs an active MATLAB license to use the source code. A "binary" stereolithography model (STL file) of the relevant geometry is also required.


Usage
--------------
Please edit `RUN_Function.m` file to enter initial conditions and change options related to various available aerothermodynamic models. The code takes in a binary STL file geometry as input, which must be specified in the `RUN_Function.m` matlab script file. Set initial conditions and type "RUN_Function" in the MATLAB command window. Make sure that you are in the same directory of the RUN file.

Aditionally, enable `opt.checks = 1` in the `RUN_Function.m` file to check the reference system, refined geometry, correct orientation of the face normals, shadowed panels, local radius quality and transitional bride function checks.

The code prints CL, CD and stagnation point heat flux as outputs for a given geometry and altitude conditions.


Under development
--------
1. In the spirit of making the software accessible to everyone, an updated version in python language is under development.
2. Multi-body object definition.
3. Re-entry trajectory analysis suite.

Disclaimer
------
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MPL license for more details.

The copyright holders are not liable for any damage(s) incurred due to improper use of `FOSTRAD`.
