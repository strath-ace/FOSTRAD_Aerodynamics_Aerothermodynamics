function [F] = ViscosityWall(Nrho, Twall)
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


Hindex = length(Twall);

%Sutherland Coefficients
%S1 is the effective temperature [K]
%O         N2          O2          Ar          N           He          H
S1 = [127       111         127         144         111         79.4        72];
S2 = [1.6934    1.4067      1.6934      21.25       1.4067      1.48438     0.63624]*1E-6;
Nrho_tot = sum(Nrho,2);

%Mean viscosity coefficient
S1mix = sum(S1.*Nrho,2)./Nrho_tot;
S2mix = sum(S2.*Nrho,2)./Nrho_tot;

VisSuWall = (S2mix .* Twall .^(3/2))./(Twall  + S1mix );

F = [VisSuWall'];