%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOSTRAD: The Free Open Source Tool for Re-entry of Asteroids and Debris %%
%%%%%%%%  Copyright (C) 2018 University of Strathclyde and Authors %%%%%%%%%%
%%%%%%%%             Aerospace Centre of Excellence,               %%%%%%%%%%
%%%%%%%%           Mechanical and Aerospace Engineering            %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% License:
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the Mozilla Public License, v. 2.0. If a copy of
%    the MPL was not distributed with this file, You can obtain one at
%    http://mozilla.org/MPL/2.0/.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    MPL license for more details.
%
%% Authors: Alessandro Falchi, Edmondo Minisci, Massimiliano Vasile,
%%          Marco Fossati, Piyush Mehta, Gianluca Benedetti, Fabio
%%          Morgado, Sai Abhishek Peddakotla
%%
%% Contact: edmondo.minisci@strath.ac.uk, massimiliano.vasile@strath.ac.uk,
%%          marco.fossati@strath.ac.uk
%
%
% Kindly acknowledge any of the relevant references where ever necessary when using FOSTRAD; 
%
%   1. Falchi, Alessandro, et al. "FOSTRAD: An advanced open source tool for re-entry analysis." 15th Reinventing Space 
%   Conference. 2017.
%   2. Mehta, Piyush, et al. "An open source hypersonic aerodynamic and aerothermodynamic modelling tool." 8th European 
%   Symposium on Aerothermodynamics for Space Vehicles. 2015.
%   3. Benedetti, Gianluca, et al. "Low-fidelity modelling for aerodynamic characteristics of re-entry objects." Stardust 
%   Final Conference. Springer, Cham, 2018.

function [F] = Adimensionals_v1(C, TFree, TWall, Velocity, d0)
% vectorialized
% Adimensionals(C0, C1, C2, C3, TFree, TWall, rho0, Velocity, d0 )
% F = [ Re, SR, Mach, Kn, Tratio]

C0 = C(:,1);
C1 = C(:,2);
C2 = C(:,3);
C3 = C(:,4);
rho0 = C(:,5);


vis0 = C0.*sqrt(TFree);  %viscosity VHS model [Pa.s]
a0 = C1.*sqrt(TFree);    %speed of sound [m/s]
MFP0 = C2./rho0;          %Mean Free Path [m]
MPS0 = C3.*sqrt(TFree); %most probable speed (for the speed ratio) [m/s]


%Adimensional numbers:
Re=rho0.*Velocity.*d0./vis0;
SR = Velocity./MPS0;
Mach=Velocity./a0;
if length(TWall) > 1
    warning('computing the Tratio for the average surface temperature')
    Tratio = mean(TWall)./TFree; 
else
    Tratio = TWall./TFree; 
end
Kn=MFP0./d0;

F = [Re, SR, Mach, Kn, Tratio];

