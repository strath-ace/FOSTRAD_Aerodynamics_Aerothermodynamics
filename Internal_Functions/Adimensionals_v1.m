function [F] = Adimensionals_v1(C, TFree, TWall, Velocity, d0)
% vectorialized
% Adimensionals(C0, C1, C2, C3, TFree, TWall, rho0, Velocity, d0 )
% F = [ Re, SR, Mach, Kn, Tratio]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Copyright (C) 2018 University of Strathclyde and Authors  %%%%%%%%
%%%%%%%          e-mail: alessandro.falchi@strath.ac.uk            %%%%%%%%
%%%%%%%                Author: Alessandro Falchi                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% License:
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%

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

