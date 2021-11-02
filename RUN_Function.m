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
%	 the MPL was not distributed with this file, You can obtain one at
% 	 http://mozilla.org/MPL/2.0/.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    MPL license for more details.
%
%% Authors: Alessandro Falchi, Edmondo Minisci, Massimiliano Vasile,
%%			Marco Fossati, Piyush Mehta, Gianluca Benedetti, Fabio
%%			Morgado, Sai Abhishek Peddakotla
%%
%% Contact: edmondo.minisci@strath.ac.uk, massimiliano.vasile@strath.ac.uk,
%%			marco.fossati@strath.ac.uk
%
%
% Kindly acknowledge any of the relevant references where ever necessary when using FOSTRAD; 
%
%	1. Falchi, Alessandro, et al. "FOSTRAD: An advanced open source tool for re-entry analysis." 15th Reinventing Space 
% 	Conference. 2017.
%	2. Mehta, Piyush, et al. "An open source hypersonic aerodynamic and aerothermodynamic modelling tool." 8th European 
%	Symposium on Aerothermodynamics for Space Vehicles. 2015.
%	3. Benedetti, Gianluca, et al. "Low-fidelity modelling for aerodynamic characteristics of re-entry objects." Stardust 
%	Final Conference. Springer, Cham, 2018.

function [FF]=RUN_Function()

format long g
%run
cDir=pwd;

addpath(genpath('./STRATH_A_main'));
addpath(genpath('./DataSets_SM'));
addpath(genpath('./External_Functions'));
addpath(genpath('./Internal_Functions'));
addpath(genpath('./Mesh_voxelisation'));
addpath(genpath('./RadiusSmoothFunc'));

%% General inputs - used by both aerodynamics and aerothermodynamics module

% Name of the .stl file to be used for the aero/aerothermodynamic characterization
% this may be used as a BINARY or ASCII format. It must be defined as a
% cell for using the multi-body definition
STLnames = {[cDir,'/IXV.00001.stl']};

% Altitude that must be simulated. It uses the NRLMSISE-00 at fixed
% latitude, longitude, day and hour of the year. The atmospheric model has
% a range from 0km to 999km, with 1km spacing, an spline interpolation is
% used to compute atmospheric parameters between each integer altitude.
% A vector of altitudes can be simulated.
altitude = 70; %[km]  min and max: 0km ~ 999km

% Velocity or speed ratio (if V < 50, it is interpreted as a molecular speed ratio)
Vinf = 6500;   %[m/s]

% Object Reference length (it can be used the equivalent diameter of the object)
% If not provided (i.e. opt.LREF = []; ), it will be automatically computed as
% an equivalent cross section diameter.
opt.LREF = [];  % [m] Reference length

% reference cross section area, used to compute aerodynamic coefficients
% If not provided (i.e. opt.SREF = []; ), it will be automatically computed.
opt.SREF = [];  % [m] Object cross section

% STL mesh scaling factor
opt.scale = 1.0;  % [1.0] by default

% Automatic computation of opt.SREF and Lref at AofA & Side-slip = 0deg
% requires a back-face culling step, useful for new geometries.
opt.autoSREFFlag = 0;

% SS and AofA must be given as single inputs. A nested for loop may be
% used to obtain multi-attitude database as shown in the example below.
AofA = -10;  % [deg]  angle of attack
SS = 0;    % [deg]  sideslip angle
Roll= 0;   % [deg]  roll angle

% initial checks, useful to check the reference system and that normals,
% geometry, local radius, etc... are being computed correctly
opt.checks = 0;

opt.fixed_COG = []; % Fixed centre of gravity coordinates [x y z];
% If empty it is automatically computed using the
% voxelator assuming the defined densities

opt.NminFaces = 100;   % minimum number of facets for the mesh refinement

opt.NmaxFaces = 500;  % maximum number of facets to be refined via the internal mesh refiner...
% If the aerothermodynamics are computed with the local radius approach
% it is suggested to test different values during the initial local radius benchmark

opt.nPix = [3000 3001]; % pixel resolution to be used for the analyses with the panel visibility algorithm
% 1500~3000 is a good compromise between accuracy and computational efficiency
% for using a convergence loop based on the visible
% area use starting and maximum pixels resolution,
% e.g.: [1500 3501];

opt.nPixStep = 250;      % Pixels increment used for the convergence loop

opt.bfcFlag = 1;        % Use the back-face culling within the facets visibility detection phase
% The panels facing backwards will be automatically
% excluded from the computation

%% Volumetric, Mass, inertial, and geometry information for a multi-body assembly
% It's used only for estimating the average properties of an object,
% determine their geometric shape evolution during a re-entry scenario.
% these properties are not used for estimating the aerodynamics and
% aerothermodynamics of the object.
% If an assembly is made of different coherent .stl files, their properties
% must be specified. In the following setup an example of an assembly made
% of 6 objects is reported.

opt.NnormGrid = 150; % voxels per unit length, used to estimate the mass of a uniform object
opt.voxFlag = 0;     % To activate the voxelization for the mass distribution
% at the first simulation it is automatically active
% in order to allow the initialization.

BodyPar.SimpleGeom = [1 1 1 1 1 1]; % Used for the normals coherence check, geometry is simple if
% ALL the normals on vertexes (or faces) and the vertexes vectors
% have a positive dot product; used for the winding direction check
% If the geometry is not simple, and the geometry has some non-coherent windings
% it should be fixed with a proper opt.MESHing software or manually.
BodyPar.rho = [2700 2700 1200 2700 1800 2700]; %[kg/m3] % Equivalent density of the simulated object
BodyPar.Q_ampl = [1 1 1 1 1 1]; % used for applying an augmented heat flux on the specified object.
BodyPar.ShellFlag = [0 0 1 1 1 0]; % if the object is a shell, it's geometric shape is not scaled with the progressive ablation
BodyPar.ThinObjFlag = [0 1 0 0 0 0]; % if the object is thin (i.e.: a flat plate) only the two major lengths
% are scaled with the progressive ablation
BodyPar.shellThickness = [0.1 0 0.005 0.005 0.005 0]; % Defined shell thickness
BodyPar.survingFaceID_sys_global = cell(1,1);         % previously used for testing a local ablation approach
BodyPar.ablFaceID_sys_global = cell(1,1);             % previously used for testing a local ablation approach



%% Module activation

%Flags used to activate the different modules.
opt.AeroFlag = 1;     % if Aeroflag = 1 the Aerodynamics module is active
opt.ThermoFlag = 1;   % if Thermoflag = 1 the Aero-Thermodynamics module is active, both aero and thermal module must be active at the same time.
opt.ablationFlag = 0;     % if opt.ablationFlag = 1 the Aero-Thermodynamics ablation module is active, the local radius is computed at each run

%% Aerodynamics module inputs

opt.cordLengthRef = [];  % Cord length used for Moment coefficient computation; if left empty it is automatically computed

% Area ratio based drag and pitching moment correction factor
opt.CF_Flag = 0; % [0] No aerodynamic correction factor applied
% [1] Parallelepiped aerodynamic correction factor applied
% [2] Cylindrical aerodynamic correction factor applied
% REF: Falchi, A., Minisci, E., Vasile, M., Rastelli, D., and Bellini, N. 2017.
% Dsmc-based correction factor for low-fidelity hypersonic aerodynamics of
% re-entering objects and space debris. In 7th European Conference for Aeronautics
% and Space Sciences

%% Aerothermodynamics module inputs

opt.RN_const_flag = 0;  % if == 1, constant reference curvature radius is used for the aerothermodynamics
opt.rN_ref = [];        % if left empty, and opt.RN_const_flag == 1, the cross section area equivalent radius is used
% if opt.RN_const_flag == 1 an arbitrary reference curvature radius may be used
opt.radiusUpdateFlag = 1; % To update the local radius at each simulation

opt.Twall = 300; % Wall temperature [K] % if empty it is automatically extracted
% from the ORION re-entry T_wall curve, loaded from dataset

opt.AccCoeff = 1; % wall molecular diffusive accommodation coefficient

% Number of smoothing point to be used by the radius smoothing algorithm, in terms of N/1000 facets
% Use the following recommended values:
% NsmoothRadius = 3; % For a hemispherical shape
% NsmoothRadius = 10; % For a geometry with rounded corners and flat faces
% NsmoothRadius = 10; % For a geometry with sharp corners and flat faces
opt.NsmoothRadius = 3;

% Limit value for weighting the flat faces. Recommended: 0.9
opt.FlatEdge = 0.9;

opt.RmaxRef = sqrt(opt.SREF/pi)*2; % maximum reference radius, in accordance to Klett heat-flux estimation on flat faces.... To be double checked.

opt.AeroThModel = 'krd';     % Select between the different models:
% SCARAB ['sc'] - > higher heat flux for (TPS/Spacecraft)
% Kemp Rose Detra ['krd'] - > super catalytic wall
% Fay Riddel ['fr'] - > fully catalytic, upper boundary
% Van Driest ['vd'] - > non catalytic, lower boundary, design for demise

%% STRATH_A_mb IXV heat transfer coefficient distribution usage example
% altitude = 60;
% AofA = -40;
% SS = 0;

opt.MESH = [];
i=1;

opt.NsmoothRadius = 10;
[strath, BodyPar, opt] = STRATH_A_mb(STLnames, altitude, Vinf, Roll, AofA, SS, BodyPar, opt);

H(:,i) = strath{1}; % Requested altitudes
Kn(:,i) = strath{2}; % Knudsen at requested altitudes
CD(:,i) = strath{3}; % Drag Coefficients
CL(:,i) = strath{4}; % Lift Coefficients
CS(:,i) = strath{5}; % Side Force Coefficients
CMx(:,i)  = strath{6}; % X-axis Moment coefficients (roll)
CMy(:,i)  = strath{7}; % Y-axis Moment coefficients (pitch)
CMz(:,i)  = strath{8}; % Z-axis moment coefficients (yaw)
if opt.ThermoFlag == 1
    StConst{1,i} = strath{9}; % St*StConst = Q [W/m2]; % Constant used to compute the heat transfer coefficient
    FVNInd = strath{10}; % Mesh matrix -
    F = FVNInd{1}; % F: Facets links indexing
    V = FVNInd{2}; % V: Vertices coordinates [x y z]
    N = FVNInd{3}; % N: facets normals
    preShadowInd = FVNInd{4}; % Visible Facets indexes before the shadowing
    Fbfc = F(preShadowInd,:); % Visible Facets indexing links
    Stmap{1,i} = strath{11};  % Visible facets Stanton Number mapping (CH)
    V(:,1) = V(:,1)-min(V(:,1)); % Nose-centered coordinates reference system
else
    F = opt.MESH{1}; % F: Facets links indexing
    V = opt.MESH{2}; % V: Vertices coordinates [x y z]
    N = opt.MESH{3}; % N: facets normals
    V(:,1) = V(:,1)-min(V(:,1)); % Nose-centered coordinates reference system
end
clc
if opt.ThermoFlag == 1
    FF=[-CL, CD, max(Stmap{i})*StConst{1,i}(1)];
else
    FF=[-CL, CD];
end
