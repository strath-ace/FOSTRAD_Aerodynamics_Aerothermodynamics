function [F] = Rotation_v1(Sideslip, AofA, roll)
%% provides the rotation matrix according to: ROT = Ry*Rz*Rx[vect]
% inputs Rotation(Sideslip,AofA,roll)
% output: [ROT] 3x3
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


    Rx = [1      0              0;
          0 cosd(roll) -sind(roll);
          0 sind(roll) cosd(roll)];          %roll
      
    Ry = [cosd(AofA) 0 sind(AofA);
          0              1         0; 
          -sind(AofA) 0 cosd(AofA)];              %pitch (Angle of Attack)
      
      
    Rz = [cosd(Sideslip) sind(Sideslip) 0;
          -sind(Sideslip) cosd(Sideslip)  0; 
          0             0        1];              %yaw (sideslip)
    ROT= Ry*Rz*Rx;

F = [ROT];
    