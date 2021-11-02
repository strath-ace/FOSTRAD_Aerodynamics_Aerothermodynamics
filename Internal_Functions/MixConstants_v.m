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

function [F] = MixConstants_v(Nrho, Tflow, VHSflag)
%% Inputs [O, N2, O2, Ar, N, He, H, Tflow,  VHSflag] in Atoms/m3
% This function is vectorized and uses the auto adaptive matrix expansion (i.e.: repmat is not required)
% for older matlab version use MixConstants
% Outputs column vector for altitudes: [C0' C1' C2' C3' RHO' P' Rmix' Cpmix' GAMMA' MMean' OMEGA' VisEC' VisSu'];
%   VHSflag; 0 for Hard Sphere viscosity || 1 for Variable HS model
%   Nrho Atomic density [O N2 O2 Ar N He H] [atoms/m3]
%   Free flow Temperature [K]
%**********************************************************************%
%the function provides as output the constants [C0 C1 C2 C3 RHO] to be used for computing the 
%values for the viscosity, speed of sound, most probable speed, and MFP as it follows:
%  RHO ----> Mass density based on the inputs [kg/m3]
%  vis = C0*sqrt(Tflow);  %viscosity VHS model [Pa.s]
%  a = C1*sqrt(Tflow);    %speed of sound [m/s]
%  MFP = C2/RHO;          %Mean Free Path [m]
%  C = C3*sqrt(Tflow); %most probable speed (for the speed ratio) [m/s]
%These values are constant for an equal relative molecular composition
%**********************************************************************%


Hindex = length(Tflow);

mN2 = 28.01340;               %molar mass of nitrogen molecule, grams/mole
mO2 = 31.99880;               %molar mass of oxigen molecule, grams/mole
mO = mO2/2;                   %molar mass of oxigen atom, grams/mole          
mN = mN2/2;                   %molar mass of Nitrogen atom, grams/mole
mAr = 39.9480;                %molar mass of Argon molecule, grams/mole
mHe = 4.0026020;              %molar mass of helium molecule, grams/mole
mH= 1.007940;                 %molar mass of Hydrogen molecule, grams/mole

%A is the matrix for Constant pressure heat capacity (Cp) coefficients for
%temperature dependance polinomial experimental law  Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species
% Cp(T) = (a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4)*R;
      %a1      a2       a3        a4        a5  
A = [3.16826710 -3.27931884E-03 6.64306396E-06 -6.12806624E-09 2.11268971E-12     %0 polinomial coefficients
     3.53100528 -1.23660987E-04 -5.02999437E-07 2.43530612E-09 -1.40881235E-12    %N2 polinomial coefficients
     3.78246636 -2.99673416E-03 9.84730200E-06 -9.68129508E-09 3.24372836E-12      %02 polinomial coefficients
     2.59316097 -1.32892944E-03 5.26503944E-06 -5.97956691E-09 2.18967862E-12     %Ar+ polinomial coefficients
     2.5        0               0               0               0                 %N 
     2.5        0               0               0               0                 %He polinomial coefficients
     2.5        0               0               0               0                 %H
     ];
 
 %Sutherland Coefficients
        %S1 is the effective temperature [K]
     %O         N2          O2          Ar          N           He          H 
S1 = [127       111         127         144         111         79.4        72];
S2 = [1.6934    1.4067      1.6934      21.25       1.4067      1.48438     0.63624]*1E-6;

 
 
%atomic densities
         %O  N2  O2  Ar  N   He  H           [kg/mol]
% Nrho =  [O N2 O2 Ar N He H];        
gamma = [5/3 1.4 1.4 5/3 5/3 5/3 5/3];
M   =   [mO  mN2 mO2 mAr mN  mHe mH]/1000;
NA=	6.0221409E+023;	%mol-1
R=	8.314472;	%J/molK
Boltzmann=	1.3806504E-23;	%J/K
Ri = R./M';

%molecularMass [kg/atom]
m = M/NA;


if VHSflag==0
    %molecular diameters HS and VHS model
    diam = [  	%O      N2          O2          Ar           N      He        H
               3E-10	3.784E-10	3.636E-10	3.659E-10    3E-10  2.33E-10  3E-10   %HS  
               %3E-10	4.17E-10	4.07E-10	4.17E-10     3E-10  2.33E-10  3E-10   %VHS
           ];
elseif VHSflag==1
    diam = [  	%O      N2          O2          Ar           NV
%                3E-10	3.784E-10	3.636E-10	3.659E-10    3E-10  2.33E-10  3E-10    %HS  
               3E-10	4.17E-10	4.07E-10	4.17E-10     3E-10  2.33E-10  3E-10    %VHS
           ];
end

    

%viscosity coefficients
         %O   N2   O2   Ar   N   He  H
omega = [0.8 0.74 0.77 0.81 0.8 0.66 0.8];
T = [ones(Hindex,1) Tflow Tflow.^2 Tflow.^3 Tflow.^4]';
temp = (A*T);
for i=1:7    
    Cp(i,1:Hindex) = Ri(i,1)*temp(i,:);
end


    N_tot = sum(Nrho,2);
    RHO = sum(Nrho.*M,2)/NA;                 %Mean mass density [kg/m3]
    DIAM = sum(diam.*Nrho,2)./N_tot;         %Mean Molecular Diameter [m]
    OMEGA = sum(omega.*Nrho,2)./N_tot;                        %Mean viscosity coefficient
    S1mix = sum(S1.*Nrho,2)./N_tot; 
    S2mix = sum(S2.*Nrho,2)./N_tot; 
    MMean = sum(M.*Nrho,2)./N_tot;           %Mean Molar Mass [kg/mol]
    mMean = sum(m.*Nrho,2)./N_tot;           %Mean Molecular Mass
    GAMMA=sum(gamma.*Nrho,2)./N_tot;       %Mean Adiabatic coefficient
    Rmix = R./MMean;
    Cpmix = sum(Cp'.*Nrho,2)./N_tot;
    P = Rmix.*Tflow.*RHO;
    
    
    
    C0 = (5*mMean/16).*sqrt(pi().*R./MMean)./(pi().*DIAM.^2);
    C1 = sqrt(GAMMA.*R./MMean);
    C2 = 2*C0./(15*sqrt(2*pi()*R./MMean)).*(5-2*OMEGA).*(7-2*OMEGA);
    C3 = sqrt(2*Rmix);

    
    VisEC = C0.*sqrt(Tflow);  %viscosity VHS model [Pa.s]
    VisSu = (S2mix .* Tflow.^(3/2))./(Tflow + S1mix);
%     VisSu = VisSu(1)*(Tflow/Tflow(1)).^(3/2).*(Tflow(1)+S1mix)./(Tflow+S1mix);
%     Vis_pl = VisSu(1).*(Tflow/Tflow(1)).^OMEGA;
    
%a = C1*sqrt(Tflow);    %speed of sound [m/s]
%MFP = C2./RHO;          %Mean Free Path [m]
%MPS = C3*sqrt(Tflow); %most probable speed (for the speed ratio) [m/s]

% SR = Vinf/MPS = Vinf/(C3*sqrt(Tflow))

%disp('   C0               C1               C2                C3        Density')
F = [C0 C1 C2 C3 RHO P Rmix Cpmix GAMMA MMean OMEGA VisEC VisSu];


