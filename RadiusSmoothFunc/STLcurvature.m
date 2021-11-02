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

function [FVNN, RadiiV_Vs_F_Fs] = STLcurvature(F, V, N, Rmin, Rmax, Smooth, ConcaveFlag)
%% summary [[F V N Nverts], {RadiiV_Vs_F_Fs}] = STLcurvature(F, V, N, Rmin, Rmax, Smooth{sphSmoothFlag, Nsmooth, 'avType', flatEdge, flatWeightFlag});
% Compute the radius on each triangular face for the given Faces and Vertices
% the radius is based on Gaussian Curvature K = k1*k2 = 1/r^2;
% the function contains also the values on the vertices. On cancave surfaces it returns a reference radius.
% Inputs: 
% - F               Nx3 connection between the points for creating the triangulated surface
% - V               Mx3 vertices points columns x,y,z respectively
% - N               Nx3 normals on faces
% - Nverts          Mx3 normals on vertices
% - Rmin =          0.05~0.2; the minimum radius expected, when the surface radius is very small R->0;
% - Rmax =          Maximum Reference Radius Value assigned to the flat vertices (i.e.: when R--> inf);
% -                 if Rmax == 0, Rmax = the average y-z plane semi-axes lenghts
% - Smooth = cell = {sphSmoothFlag, Nsmooth, avType} where:
% - sphSmoothFlag = if [==1] Spherical Searchable Volume Smoothing enabled: find the indexes within a sphere
%                   built on each vertex of defined radius based on the number of smoothing points required that
%                   the indexes are sorted according to the distance, into a descending order.
%                   dramatically increases the computational times for huge meshes; if the number of 
%                   vertices is above 30k the parallel processing is enabled);
%                   if [==0] the averaged gaussian radius non-smoothed is given as output.
% - Nsmooth = Nsmooth is the number of averaged points. Increase or decrease the number of averaged points
%                   for increasing or decreasing the smoothness of the radius, the maximum Nsmooth is given by length(V)
% - avType  = char 'e' for an exponential moving average; 's' for a simple average; 
% - flatEdge= the algorithm detects the number of flatVertices, according to the ratio of flatVertices/Nsmooth within the 
%                       selected smoothing points, it assign the reference flat radius to the vertex.
%                       A vertex is considered flat when flatVertices/Nsmooth >= flatEdge.
%                   *** if flatEdge == 1, only the Vertices whose ALL nearest vertexes are "flat" will be considered flat
%                       all the remaining vertices (the ones with flatVertices/Nsmooth < 1) will be averaged.
%                   *** if flatEdge == 0, all vertices are given the flatRadiusReferenceValue, and the spherical smoothing will be
%                       DISABLED.
%                   *** if 0 < flatEdge < 1 it averages the vertices only when their flatVertices/Nsmooth < flatEdge
%                       in general, a higher flatEdge will provide a smoother result, decreasing the number of
%                       assigned flatReferenceRadius.
% - flatWeightFlag: Enables the weighting of flatVertices Reference Radius value according to the flatVertices/Nsmooth ratio on each point.
%                   *** if == 1, linear weighting of the Reference Radius 
%                   *** if == 2, sigmoid weighting of the Reference Radius (smaller edge angles, better smoothing toward the flat surfaces)
% - ConcaveFlag     if == 1, assign a flat ref radius to the concave-corners, otherwise, it computes the radius (useful for reppresenting the shock impingement e.g.: on IXV flaps)
%                   by default it is turned on.
%
% Outputs:
% - FVN {cell} 1x3 containing the reordered and cleaned triangulation {F, V, N}
% - RadiiV_Vs_F_Fs {cell} 1x4 containing: {radiiOnVerts, radiiOnVertsSmooth, radiiOnFaces, radiiOnFacesSmooth}:
% - radiiOnVerts-Smooth = matrix Mx1, where M is the number of "reordered and cleaned triagulation" Vertices
% - radiiOnFaces-Smooth = matrix Nx1, where N is the number of "reordered and cleaned triagulation" Faces
%
% Example for calling the function: 
% [V,F,N] = import_stl_fast('ASCIIstlName.stl');
% [F,V,N] = stlread('Binary_stlName.stl');
% Smooth = {1, 20, 'e', 0.9, 2};
% STLcurvature(F, V, 0.1, 0, Smooth);
%
% Note: The function employs the corner voronoi weights for averaging from vertices to face properties in accordance with the
% triangular area computed with Heron's formula. The Local curvature radius is basedo on a modified version of the algorithm
% defined and Implemented according to "Estimating Curvatures and Their Derivatives on Triangle Meshes" by Szymon Rusinkiewicz (2004)"



switch iscell(Smooth)
    case 1
        if Smooth{1} ~= 0;
            sphSmoothFlag = Smooth{1};
            Nsmooth = Smooth{2};
            if Nsmooth < 1
                Nsmooth = 1;
            end
            avType = Smooth{3};
            flatEdge = Smooth{4};
            flatWeightFlag = Smooth{5};
            if strcmp(avType, 'e') ~= 1 && strcmp(avType, 's') ~= 1
                error('avType not valid. Valid avType values: ''e'' exponential moving average and ''s'' simple average');
            end
        else
            sphSmoothFlag = 0;
        end
    otherwise
        sphSmoothFlag = 0;
end

if Rmax == 0
    disp('************************************** ATTENTION ***************************************************')
    disp('When the surface is flat a Reference [MaxRefRadius] is given as the mean Z-Y semi-axis length, it should be')
    disp('computed in accordance with a reference cross section area, depending on the "view" direction')
    disp('****************************************************************************************************')
end
% STL READ example: binary
%  [F,V,N] = stlread('IXV_original.STL')
% STL READ example: ascii
%  [V,F,N] = import_stl_fast('IXV_original.STL')


%%  Re ordering and cleaning mesh
% [V, F] = reOrderMesh(V, F, Nround);
Nvert = length(V);

% trisurf(F,V(:,1),V(:,2),V(:,3),'FaceColor','w')

if Rmax == 0
    MaxRefRadius = round(mean(max(V(:,2:3))-min(V(:,2:3)))/2,7);
else
    MaxRefRadius = Rmax;
end


% computing normals and curvatures
[CurvUxVy, ~ , Avertex, Acorner, Nverts, avEdge ] = GetCurvatures( F, V, N, 0);
MaxCurv = max(CurvUxVy,[],2);
MaxCurv(MaxCurv < 0) = 0;
MinCurv = min(CurvUxVy,[],2);
MinCurv(MinCurv > 0) = 0;


% negInd = find( [abs(MaxCurv) < abs(MinCurv)] == 1 );
% posInd = setdiff([1:Nvert]', negInd );

Sign = ones(length(CurvUxVy),1);
% Sign(CurvUxVy(:,1) < 0) = -1;
Sign((abs(MaxCurv) < abs(MinCurv)) == 1) = -1;

% Sign(CurvUxVy(:,1) < 0) = -1;
% Sign(CurvUxVy(:,2) < 0) = -1;
% NegEdge = abs(mean(CurvUxVy(CurvUxVy<0)));
% Sign(abs(CurvUxVy(:,1)) < NegEdge/10) = 1;
% Sign(abs(CurvUxVy(:,2)) < NegEdge/10) = 1;

% Test(CurvUxVy(:,1) < 1/MaxRefRadius & CurvUxVy(:,2) < 1/MaxRefRadius) = 1;
CurvUxVy = abs(CurvUxVy);
radiiOnVerts = sqrt(1./(CurvUxVy(:,1).*CurvUxVy(:,2)));


% when only one Curvature is --> 0 (flat vertex), select the non-zero curv instead of the mean curvature 
radiiOnVerts(CurvUxVy(:,1) > 1/MaxRefRadius & CurvUxVy(:,2)< 1/MaxRefRadius,1) = 1./CurvUxVy(CurvUxVy(:,1) > 1/MaxRefRadius & CurvUxVy(:,2)< 1/MaxRefRadius,1,1);
radiiOnVerts(CurvUxVy(:,2) > 1/MaxRefRadius & CurvUxVy(:,1)< 1/MaxRefRadius,1) = 1./CurvUxVy(CurvUxVy(:,2) > 1/MaxRefRadius & CurvUxVy(:,1)< 1/MaxRefRadius,2);

% outliers = find(radiiOnVerts>=MaxRefRadius);
% validInd = [1:1:length(V)]';
% % validInd(outliers) = 0;
% validInd1 = [1:1:length(V)]';
% radii0 = find(radiiOnVerts < Rmin);
% validInd1(radii0) = 0;

% AvValidRadii = mean(radiiOnVerts(validInd~=0 & validInd1~=0 & isnan(radiiOnVerts)==0));
radiiOnVerts(radiiOnVerts>=MaxRefRadius) = MaxRefRadius;
radiiOnVerts(radiiOnVerts < Rmin) = Rmin;
radiiOnVerts(isnan(radiiOnVerts)) = Rmin;

if ~exist('ConcaveFlag', 'var') 
    radiiOnVerts(Sign(:,1) < 0 ) = MaxRefRadius;
elseif ConcaveFlag == 1
    radiiOnVerts(Sign(:,1) < 0 ) = MaxRefRadius;
end




%% Searchable Sphere Algorithm

if sphSmoothFlag == 1;
    tstart = tic;
    
    radiiOnVertsSmooth = zeros(Nvert,1);
    radiiOnVertsS = single(radiiOnVerts);
%     F = single(F);
    SearchRadius = single(avEdge);
    ParFlag = Nvert/50000;
    
    if ParFlag >= 1
        disp('Starting the Parallel Local Radius Smoothing')
        parfor i = 1:Nvert
            %         InnerPointsInd(i,1) = {searchableRadius(V(i,:), V, SearchRadius, Nsmooth)};
            InnerPointsInd = {searchableRadius(V(i,:), V, SearchRadius, Nsmooth)};
            radiiOnVertsSmooth(i,1) = sphVolSmoothing(InnerPointsInd, {radiiOnVertsS}, avType, Nsmooth, MaxRefRadius, flatEdge, flatWeightFlag);
        end
    else
        disp('Starting the Local Radius Smoothing')
        for i = 1:Nvert
            %         InnerPointsInd(i,1) = {searchableRadius(V(i,:), V, SearchRadius, Nsmooth)};
            InnerPointsInd = {searchableRadius(V(i,:), V, SearchRadius, Nsmooth)};
            radiiOnVertsSmooth(i,1) = sphVolSmoothing(InnerPointsInd, {radiiOnVertsS}, avType, Nsmooth, MaxRefRadius, flatEdge, flatWeightFlag);    
        end
    end
    disp(['Surface Local Radius Smoothed in: ',num2str(toc(tstart)),'s'])
end

%%
%calculate voronoi weights based on heron's formula area.
wfp(:,1)=Acorner(:,1)./Avertex(F(:,1));
wfp(:,2)=Acorner(:,2)./Avertex(F(:,2));
wfp(:,3)=Acorner(:,3)./Avertex(F(:,3));
radiiOnFaces= sum(radiiOnVerts(F).*wfp,2)./sum(wfp,2);

% % Voronoi weights baricenter
% faceBaricenter(:,1) = sum(reshape(V(F(:,1:3),1),[],3).*wfp,2)./sum(wfp,2);
% faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3).*wfp,2)./sum(wfp,2);
% faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3).*wfp,2)./sum(wfp,2);


if sphSmoothFlag == 1
    radiiOnFacesSmooth= sum(radiiOnVertsSmooth(F).*wfp,2)./sum(wfp,2);
    RadiiV_Vs_F_Fs = [{radiiOnVerts}, {radiiOnVertsSmooth}, {radiiOnFaces}, {radiiOnFacesSmooth}];
else
    RadiiV_Vs_F_Fs = [{radiiOnVerts}, {radiiOnVertsSmooth}, {radiiOnFaces}, {radiiOnFaces}];
end

% outputs
FVNN = [{F}, {V}, {N}, {Nverts}];

return
%% plotting

% % normal baricenter computation
faceBaricenter(:,1) = sum(reshape(V(F(:,1:3),1),[],3),2)/3;
faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3),2)/3;
faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3),2)/3;
% 
% maxRF = max(radiiOnFaces);
% minRF = min(radiiOnFaces);
% % 
% % Original on Face
% figure('Position', [100, 100, 1049, 895]);
% trisurf(F,V(:,1),V(:,2),V(:,3),radiiOnFaces,'EdgeAlpha',0.15);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Average Radius [m]')
% set(gca, 'CLim', [min(radiiOnFaces), max(radiiOnFaces)]);
% view([-45 5])
% colorbar
% 
% % Smoothed
% if sphSmoothFlag == 1
%     figure('Position', [100, 100, 1049, 895]);
%     trisurf(F,V(:,1),V(:,2),V(:,3),radiiOnFacesSmooth,'EdgeAlpha',0.15);
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title('IXV Smoothed Local Radius [m]')
%     set(gca, 'CLim', [min(radiiOnFacesSmooth), max(radiiOnFacesSmooth)]);
% %     set(gca, 'CLim', [0.5, 2]);
%     view([-45 5])
%     colorbar
% end
% 
% 
% MaxC = max(abs((radiiOnFacesSmooth-radiiOnFaces)./radiiOnFacesSmooth*100));
% % % plotting the difference between smoothing and non smooting
% figure('Position', [100, 100, 1049, 895]);
% trisurf(F,V(:,1),V(:,2),V(:,3),(radiiOnFacesSmooth-radiiOnFaces)./radiiOnFacesSmooth*100,'EdgeAlpha',0.15);
% % trisurf(F,V(:,1),V(:,2),V(:,3),(radiiOnFacesSmooth2-radiiOnFacesSmooth)./radiiOnFacesSmooth2*100,'EdgeAlpha',0.15);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % title('IXV Smoothed Local Radius [m]')
% view([-45 5])
% colorbar;
% c= colorbar;
% set(gca, 'CLim', [-MaxC MaxC]);
% % set(gca, 'CLim', [-0.01 0.01]);
% ylabel(c,'Local Radius [m] ','FontSize',16)
% title(['Sigmoid Smoothed vs Non-smoothed, Nsmooth = ',num2str(Nsmooth)],'FontSize',18)
% set(gca, 'FontSize', 14, 'FontName','Sans')
% 
% 
% figure()
% view([-37.5 30])
% grid on
% hold on
% trisurf(F,V(:,1),V(:,2),V(:,3),radiiOnFaces,'EdgeAlpha',0.15)
% % fq1 = quiver3(faceBaricenter(:,1),faceBaricenter(:,2),faceBaricenter(:,3), -N(:,1), -N(:,2), -N(:,3),3,'r');
% % set(fq1,'linewidth',1.5);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Normals on faces')
% set(gca, 'CLim', [minRF, maxRF]);
% colorbar
% 
figure()
view([-37.5 30])
grid on
hold on
trisurf(F,V(:,1),V(:,2),V(:,3),'FaceColor','w') 
scatter3(V(:,1),V(:,2),V(:,3),30,radiiOnVerts,'o','LineWidth',20)
 fq2 = quiver3(V(:,1),V(:,2),V(:,3), Nverts(:,1),Nverts(:,2),Nverts(:,3),2,'r');
% set(fq2,'linewidth',1.5);
xlabel('x')
ylabel('y')
zlabel('z')
title('Normals and Radius on Vertices')
set(gca, 'CLim', [minRF, maxRF]);
colorbar
end
