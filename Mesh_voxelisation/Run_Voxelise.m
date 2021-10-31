function [V_vox, voxBar, I_vox, dVol, mass] = Run_Voxelise(BodyPar, F, V, NnormGrid, rho, ShellFlag, shellThickness, bodyID )
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
% computes the voxelization of an arbitrary mesh defined by the faces [F] and vertexes [V]
% returns also the intertia tensor, provided a uniform density rho

global bbox innerFlag

% checking if it's running a Multi body simulation
if ~isempty('BodyPar')
    mbFlag = 1;
else
    mbFlag = 0;
end

if ~exist('bodyID','var')
    bodyID = 1;
end

tstart = tic;
if ~exist('rho','var')
    disp('Computing the normalized Inertia Tensors: rho = 1kg/m3')
    rho = 1;
end

if ~exist('ShellFlag','var')
    ShellFlag = 0;    
end

if ShellFlag ==1
    if ~exist('shellThickness','var')
        shellThickness = input('Please input The shell thickness in [m]: >> ');
    end
end

% % substituted with the voxel-based baricenter
% faceBaricenter = (V(F(:,1),1)+V(F(:,2),1)+V(F(:,3),1))/3;
% faceBaricenter(:,2) = (V(F(:,1),2)+V(F(:,2),2)+V(F(:,3),2))/3;
% faceBaricenter(:,3) = (V(F(:,1),3)+V(F(:,2),3)+V(F(:,3),3))/3;
% A = Area(F,V);
% CoG = sum(faceBaricenter.*repmat(A,1,3))/sum(A);
% V(:,1) = V(:,1)-CoG(1);
% V(:,2) = V(:,2)-CoG(2);
% V(:,3) = V(:,3)-CoG(3);



% creating the mesh structure
if ShellFlag == 0
    mesh.faces = F;
    mesh.vertices = V;

    % finding the bounding box
    bbox = [min(V); max(V)];
    bLim = [bbox(2,1)-bbox(1,1), bbox(2,3)-bbox(1,3), bbox(2,3)-bbox(1,3)];
    Ngrid = ceil(NnormGrid.*bLim./max(bLim));
    dxyz = bLim./Ngrid;
    dVol = dxyz(1)*dxyz(2)*dxyz(3);

    % running the voxelization
    [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(Ngrid(1), Ngrid(2), Ngrid(3), mesh, ShellFlag);

    [X, Y, Z] = ndgrid(gridCOx,gridCOy,gridCOz);
%     X1 = X(gridOUTPUT == 0);
%     Y1 = Y(gridOUTPUT == 0);
%     Z1 = Z(gridOUTPUT == 0); 
%     V_vox1 = [X1, Y1, Z1];    
    X = X(gridOUTPUT == 1);
    Y = Y(gridOUTPUT == 1);
    Z = Z(gridOUTPUT == 1);    
    V_vox = [X, Y, Z];
    
else        
    bbox = [];
    innerFlag = [];   
    % outer mesh
    faceBaricenter = sum(reshape(V(F(:,1:3),1),[],3),2)/3;
    faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3),2)/3;
    faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3),2)/3;
    Afaces = Area(F,V);
    COG = sum(faceBaricenter.*repmat(Afaces,1,3))./sum(Afaces);
    V = V - repmat(COG,length(V),1);
    mesh.faces = F;
    mesh.vertices = V;    
    
    % finding the bounding box
    bbox = [min(V); max(V)];
    bLim = [bbox(2,1)-bbox(1,1), bbox(2,3)-bbox(1,3), bbox(2,3)-bbox(1,3)];
    Ngrid = ceil(NnormGrid.*bLim./max(bLim));
    dxyz = bLim./Ngrid;
    dVol = dxyz(1)*dxyz(2)*dxyz(3);      
    if shellThickness < dxyz(1) && ~(isfield(BodyPar,'V_InInnerShell') && length(BodyPar.V_InInnerShell) >= bodyID)
        error('The Shell is too Thin to be computed with the voxelator, increase the refinement of the grid')
    end
    
    % creating the inner shell
    if mbFlag == 1 && ~isfield(BodyPar,'V_InInnerShell')
        V_s1 = V(:,1).*(1-shellThickness*2/bLim(1));
        V_s1(:,2) = V(:,2).*(1-shellThickness*2/bLim(2));
        V_s1(:,3) = V(:,3).*(1-shellThickness*2/bLim(3));
        % outer mesh
        mesh_s1.faces = F;
        mesh_s1.vertices = V_s1;  
        BodyPar.V_InInnerShell{bodyID} = V_s1;
    elseif mbFlag == 1 && isfield(BodyPar,'V_InInnerShell') && length(BodyPar.V_InInnerShell) < bodyID
        V_s1 = V(:,1).*(1-shellThickness*2/bLim(1));
        V_s1(:,2) = V(:,2).*(1-shellThickness*2/bLim(2));
        V_s1(:,3) = V(:,3).*(1-shellThickness*2/bLim(3));
        % outer mesh
        mesh_s1.faces = F;
        mesh_s1.vertices = V_s1;  
        BodyPar.V_InInnerShell{bodyID} = V_s1;
    elseif mbFlag == 1 && isfield(BodyPar,'V_InInnerShell') && length(BodyPar.V_InInnerShell) >= bodyID
        mesh_s1.faces = F;
        mesh_s1.vertices = BodyPar.V_InInnerShell{bodyID};  
    elseif mbFlag == 0
        error('This version is supposed to be used with the Multi-body simulation: try the Run_Voxelise_old.m')
    end
    
    % running the outer shell voxelization
    [gridOUTPUT] = VOXELISE(Ngrid(1), Ngrid(2), Ngrid(3), mesh, ShellFlag);
   
    % running the inner shell voxelization
    [gridOUTPUT_s1,gridCOx,gridCOy,gridCOz] = VOXELISE(Ngrid(1), Ngrid(2), Ngrid(3), mesh_s1, ShellFlag);
    [X, Y, Z] = ndgrid(gridCOx,gridCOy,gridCOz);
    X = X(gridOUTPUT == 1 & gridOUTPUT_s1 == 0);
    Y = Y(gridOUTPUT == 1 & gridOUTPUT_s1 == 0);
    Z = Z(gridOUTPUT == 1 & gridOUTPUT_s1 == 0);
    V_vox = [X, Y, Z]+repmat(COG,length(X),1);   
end

% for shells the baricenter is more accurate, given that they are simple geometries
if ShellFlag == 1
    voxBar = COG;
else
    voxBar = mean(V_vox);
end
V_vox = V_vox - repmat(voxBar,length(V_vox),1);
mass = dVol*rho*length(V_vox);

Ixx = sum((V_vox(:,2).^2 + V_vox(:,3).^2)*dVol*rho);
Iyy = sum((V_vox(:,1).^2 + V_vox(:,3).^2)*dVol*rho);
Izz = sum((V_vox(:,1).^2 + V_vox(:,2).^2)*dVol*rho);
Ixy = 0; % -sum((V_vox(:,1).*V_vox(:,2))*dVol*rho);
Ixz = 0; % -sum((V_vox(:,1).*V_vox(:,3))*dVol*rho);
Iyz = 0; % -sum((V_vox(:,2).*V_vox(:,3))*dVol*rho);

I_vox = [Ixx Ixy Ixz
        Ixy Iyy Iyz
        Ixz Iyz Izz];

disp(['Voxelization computed in: ', num2str(toc(tstart)),'s'])

end


