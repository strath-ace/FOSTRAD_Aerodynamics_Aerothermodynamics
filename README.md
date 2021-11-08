# FOSTRAD

(The Free Open Source Tool for Re-entry of Asteroids and Debris)

`FOSTRAD` is a simulation suite that allows for the estimation of aerodynamics and aerothermodynamics of an object experiencing flow-fields during atmospheric re-entry scenarios. It has better prediction capability than object-oriented re-entry codes as it employs local panel formulation to compute aerodynamic and aerothermodynamic coefficients. A


Requirements
------

`FOSTRAD` is currently written in `MATLAB` and hence needs an active MATLAB license to use the source code. A "binary" STL file of the relevant geometry is also required.


Usage
--------------
The code takes in a binary STL file geometry as input, which must be specified in the `RUN_Function.m` matlab script file.

Aerothermodynamic models:


Under development
--------
1. In the spirit of making the software accessible to everyone, an updated version in python language is under development.
2. Multi-body object definition.
3. Re-entry trajectory analysis suite.

Disclaimer
------
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MPL license for more details.

The copyright holders are not liable for any damage(s) incurred due to improper use of `FOSTRAD`.
