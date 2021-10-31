function [strath, BodyPar, opt] = STRATH_A_mb(STLnames, altitude, Vinf, Roll, AofA, SS, BodyPar, opt) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STRATH-A: Strathclyde Tool for Re-entry and Aero-THermal Analysis %%%%
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
%% Inputs have been described in the < Controller_STRATH_mb.m > file

%% Outpus:
% [strath, BodyPar, opt] = STRATH_A_mb(STLnames, altitude, Vinf, Roll, AofA, SS, BodyPar, opt);  
% 
% H(:,i) = strath{1}; % Requested altitudes
% Kn(:,i) = strath{2}; % Knudsen at requested altitudes
% CD(:,i) = strath{3}; % Drag Coefficients
% CL(:,i) = strath{4}; % Lift Coefficients
% CS(:,i) = strath{5}; % Side Force Coefficients
% CMx(:,i)  = strath{6}; % X-axis Moment coefficients (roll)
% CMy(:,i)  = strath{7}; % Y-axis Moment coefficients (pitch)
% CMz(:,i)  = strath{8}; % Z-axis moment coefficients (yaw)
% if opt.ThermoFlag == 1
%     StConst{1,i} = strath{9}; % St*StConst = Q [W/m2]; % Constant used to compute the heat transfer coefficient
%     FVNInd = strath{10}; % Mesh matrix - 
%     F = FVNInd{1}; % F: Facets links indexing
%     V = FVNInd{2}; % V: Vertices coordinates [x y z]
%     N = FVNInd{3}; % N: facets normals
%     preShadowInd = FVNInd{4}; % Visible Facets indexes before the shadowing
%     Fbfc = F(preShadowInd,:); % Visible Facets indexing links
%     Stmap{1,i} = strath{11};  % Visible facets Stanton Number mapping (CH)
%     V(:,1) = V(:,1)-min(V(:,1)); % Nose-centered coordinates reference system
% else
%     F = opt.MESH{1};
%     V = opt.MESH{2};
%     N = opt.MESH{3};
%     V(:,1) = V(:,1)-min(V(:,1));
% end

%%


% The shadowing refinement functionality has not been tested within this
% version, and it should be used only for specific attitude/altitude based
% simulations.
shadowRefFlag = 0;


%% %%%%%%%***************  Additional inputs  ********************%%%%%%%%%
tstart = tic;


% ******************* Additional User Inputs ********************
% hmin=min(altitude); %km              %minimum altitude present in the database: 0km
% hmax=max(altitude); %km              %maximum altitude present in the database: 999km
if altitude > 999
    warning('Maximum Altitude Exceeded, H = 998km')
    altitude(altitude > 999) = 998;
end
if altitude < 1
    warning('Minimum Altitude Exceeded, H = 1km')
    altitude(altitude < 1) = 1;
end

Hcont = 0;             %Fixed reference continuum altitude
Hfmf = 500;            %Fixed reference Free molecular flow altitude
limKn_inf= 1.0E-4;     %fixed reference Kn continuum limit % affects the transition between transitional and continuum (Calibrated for the Erf Bridging)
limKn_sup= 100;        %fixed reference Kn FMF limit
limKn_inf_heat = 1e-3; %reference for the heat transfer based on the local radius approach


M2M_RR = 10; % Max to Min local radius ratio

% transitional regime calibration on the error function
CF_ratiolow = 1e-6;   % calibration factor toward the continuum regime
CF_ratiohigh = 0.1508;   % calibration factor toward the free molecular regime obtained via the calibration on 
                         % the data of the Mars Pathfinder, Mars Micro Probe, Orion CEV, Apollo CM, STS Orbiter and Sphere
% to use the normal Erf just use a very low CF = 1e-6;

% load the stagnation point heating bridging function based on Mars Pathfinder, Mars Micro Probe, Orion CEV, and Sphere:
% load('PchipHeatingBridge.mat');  % loads the BridgeCoeffPchip coefficients % obsolete 
% load('Bridging_Rn_Kn_SM.mat');
load('Bridging_Rn_Kn_SM_v2.mat'); %% Contains the additional blunted cone 5deg ref. moss 1987

rhosl = 1.225; % density at sea level [kg/m^3]

% AofA = 0;

load('OrionDatabase')  % used for the Wall Temperature estimation when the temperature is not provided
if length(altitude) > 1
    Tw = pchip(OrionDatabase(:,1), OrionDatabase(:,2), median(altitude));
else
    if altitude < 23
        Tw = 300;
    else
        Tw = pchip(OrionDatabase(:,1), OrionDatabase(:,2), altitude);
    end
end

if ~isempty(opt.Twall)
    Tw = opt.Twall;
end

% % Tw = 1000;                %fixed wall temperature

% Atmospheric model and GSI constant inputs
% opt.AccCoeff = 1;
SN = opt.AccCoeff;          %Normal Momentum Accommodation coefficient
ST = opt.AccCoeff;          %Tangential Momentum Accommodation Coefficient
AT = ST*(2-ST);         %Tangental Energy Accomodation Coefficient
AN = 1 - (1 - SN)^2;    %Normal Energy Accommodation Coefficient
AC = 0.5 * (AN + AT);   %Overall Energy Accommodation Coefficient
VHSflag = 0;            % Variable Hard Sphere Model for MFP computation VHSflag == 1, Hard Sphere model with VHSflag = 0;
if ~exist('opt.AeroThModel','var') || isempty(opt.AeroThModel)
    opt.AeroThModel = 'sc';     % Select between the different models: 
                        % SCARAB ['sc'] - > higher heat flux for (TPS/Spacecraft)
                        % Kemp Rose Detra ['krd'] - > super catalytic catalytic
                        % Fay Riddel ['fr'] - > fully catalytic, upper boundary
                        % Van Driest ['sc'] - > non catalytic, lower boundary, design for demise
end

% Shadowing algorithm convergence loop inputs
Nref = 1.5; % number of times the shadow needs to be refined
if ~isfield(opt,'nPixStep') || isempty(opt.nPixStep)
    opt.nPixStep = 250;
end
n_pix_temp = opt.nPix;
if length(n_pix_temp) == 2
    maxPixRes = n_pix_temp(2);
    n_pix_temp = n_pix_temp(1);
else
    maxPixRes = n_pix_temp+opt.nPixStep*4+1;
end
AreaConvError = 0.0025;


% NRLMSISE-00 DATABASE condition for:
% Year= 2000, Month=  1, Day=  1, Hour= 1.50,
% Time_type = Universal
% Coordinate_type = Geographic
% Latitude=   55.00deg, Longitude=   45.00deg,
% Database format [1000x9]
% Selected output parameters:
% 1  Height, kmmaxFaces
% 2  O, cm-3
% 3  N2, cm-3
% 4  O2, cm-3
% 5  Temperature_neutral, K
% 6  He, cm-3
% 7  Ar, cm-3
% 8  H, cm-3
% 9  N, cm-3
load('NRLMSISE00')
%A database with different condition may be generated and download from:
%     http://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php


Nbodies = length(STLnames);
if length(opt.NmaxFaces) < Nbodies
    opt.NmaxFaces = opt.NmaxFaces(1).*ones(1,Nbodies);
end

if any([length(BodyPar.rho),  length(BodyPar.ShellFlag), length(BodyPar.shellThickness)] < Nbodies)
    error('Check the materials, and shell properties for each body part')
end

if opt.checks == 1
    opt.checks = input('Do you want to continue the initial checks? Reply [0 1] >> ');
end

% Loading triangular opt.MESH
if ~isfield(opt,'MESH') || isempty(opt.MESH)
    for i = 1:Nbodies
        if isbinary(STLnames{i},100)
            [BodyPar.F{i}, BodyPar.V_sys{i}, ~] = stlread(STLnames{i});      % Binary stl
        else
            [BodyPar.V_sys{i}, BodyPar.F{i}, ~] = import_stl_fast(STLnames{i},1);       % ASCII stl
        end
        
        if isfield(opt,'scale') && ~isempty(opt.scale)
            BodyPar.V_sys{i} = BodyPar.V_sys{i}.*opt.scale;
        end
        
        Afaces = Area(BodyPar.F{i},BodyPar.V_sys{i});        
        Atot(i) = sum(Afaces);
        Amean(i,1) = Atot(i)./opt.NmaxFaces(i);
%         if i > 1
%             F = [F; F_temp];
%             V = [V; V_temp];
%         else
%             F = F_temp;
%             V = V_temp;
%         end
%         boundBox(i,:) = max(BodyPar.V_sys{i})-min(BodyPar.V_sys{i});
%         boundBoxVol(i,1) = boundBox(i,1)*boundBox(i,2)*boundBox(i,3);
%         opt.cordLengthRef(i,1) = max(BodyPar.V_sys{i}(:,1))-min(BodyPar.V_sys{i}(:,1));  
        
        if opt.checks == 1 && isempty(opt.MESH)
            fig = figure();
            fig.Position = [100, 25, 950, 700];
            trisurf(BodyPar.F{i},BodyPar.V_sys{i}(:,1),BodyPar.V_sys{i}(:,2),BodyPar.V_sys{i}(:,3),'edgealpha',0.25,'FaceColor',...
                [0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
            disp(['Load Triangular opt.MESH Check']);
            keyboard()
            close
        end               
    end    
    

    % trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
    % trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
    % stdPlot(Labels,Title)
    % AreaTest = sum(Area(F,V))
    % Aratio = opt.SREF/AreaTest
    % [l1 l2 l3]
    


    for i = 1:Nbodies
        % opt.LREF = max(max(V)-min(V));  
%         Afaces = Area(BodyPar.F{i},BodyPar.V_sys{i});
        % Astd = std(Afaces)./length(F);
        opt.NmaxFaces = Atot/max(Amean);
        NrefIter = ceil(log2(opt.NmaxFaces(i)/length(BodyPar.F{i})));
        Amax = sum(Afaces)./(length(BodyPar.F{i})*2^(NrefIter));
        %Amax=Amax/1000;
        % Amax = 0.040;
        Lmax = sqrt(max(Amax)/pi);
        STopt.LREFFlag = 1;       

        if STopt.LREFFlag == 1
            [BodyPar.F{i}, BodyPar.V_sys{i}] = refineSTL_tri_v2_2(BodyPar.F{i}, BodyPar.V_sys{i}, 'area',[Amax],max([opt.NminFaces opt.NmaxFaces(i)]),300,{'on', (1:1:size(BodyPar.F{i},1))});
        end         
        
        % cleaning the opt.MESH
        [BodyPar.F{i}, BodyPar.V_sys{i}] = cleanVerts(BodyPar.F{i}, BodyPar.V_sys{i});

        % computing the centre of gravity and the mass for each body part
        [BodyPar.V_vox_cog{i}, BodyPar.COG(i,:), BodyPar.I{i}, BodyPar.dVol(i,1), BodyPar.mass(i,1)] = ...
                                        Run_Voxelise(BodyPar, BodyPar.F{i}, BodyPar.V_sys{i}, opt.NnormGrid, BodyPar.rho(i), BodyPar.ShellFlag(i), BodyPar.shellThickness(i), i);
        BodyPar.mass_initial(i,1) = BodyPar.mass(i,1);
        % By default the component voxelization is based on the body CofG-fixed frame
        BodyPar.V_cog{i} = BodyPar.V_sys{i} - repmat(BodyPar.COG(i,:), length(BodyPar.V_sys{i}),1);
        
        % Computing equivalent thickness:
        if BodyPar.ShellFlag(i) == 1
            BodyPar.VirtualThickness(i) = BodyPar.shellThickness(i);
        elseif BodyPar.ThinObjFlag(i) == 1
            BodyPar.VirtualThickness(i) = min(max(BodyPar.V_sys{i})-min(BodyPar.V_sys{i}))/2;
        else % equivalent radius
            BodyPar.VirtualThickness(i) = (BodyPar.dVol(i,1)*length(BodyPar.V_vox_cog{i})*3/(4*pi)).^(1/3);
        end
        
        % building the assembly
        if i == 1
            F = BodyPar.F{i};
            V = BodyPar.V_sys{i};
            BodyPar.rhoFaces = BodyPar.rho(i).*ones(length(BodyPar.F{i}),1);
            BodyPar.rhoVerts = BodyPar.rho(i).*ones(length(BodyPar.V_sys{i}),1);
            BodyPar.LimtsF(1,1:2) = [1 length(BodyPar.F{i})];
            BodyPar.LimtsV(1,1:2) = [1 length(BodyPar.V_sys{i})];                      
        else
            F = [F; BodyPar.F{i}+BodyPar.LimtsV(i-1,2)];
            V = [V; BodyPar.V_sys{i}];
            BodyPar.rhoFaces = [BodyPar.rhoFaces; BodyPar.rho(i).*ones(length(BodyPar.F{i}),1)];
            BodyPar.rhoVerts = [BodyPar.rhoVerts; BodyPar.rho(i).*ones(length(BodyPar.V_sys{i}),1)];
            BodyPar.LimtsF(i,1:2) = [BodyPar.LimtsF(i-1,2)+1 BodyPar.LimtsF(i-1,2)+length(BodyPar.F{i})];
            BodyPar.LimtsV(i,1:2) = [BodyPar.LimtsV(i-1,2)+1 BodyPar.LimtsV(i-1,2)+length(BodyPar.V_sys{i})];  
        end  
        
    end
    Afaces = Area(F,V);
    % storing the initial mass
    BodyPar.InMass = sum(BodyPar.mass);
    
    if opt.checks == 1
        fig = figure();
        fig.Position = [100, 25, 950, 700];
        trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
        plot3(BodyPar.COG(:,1),BodyPar.COG(:,2),BodyPar.COG(:,3),'or','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','r')
        disp('Refined Triangular opt.MESH Check');
        keyboard
        close
    end
    
    BodyPar.Isys = zeros(3,3);
    BodyPar.COGsys = sum(BodyPar.COG.*repmat(BodyPar.mass,1,3),1)./sum(BodyPar.mass);
    BodyPar.COGsys_hist(1,:)  = BodyPar.COGsys;
    BodyPar.COG = BodyPar.COG-repmat(BodyPar.COGsys, Nbodies,1);
    BodyPar.COG_hist{1,1}  = BodyPar.COG;
    for i = 1:Nbodies
        BodyPar.V_sys{i} = BodyPar.V_sys{i}-repmat(BodyPar.COGsys, length(BodyPar.V_sys{i}),1);        
        BodyPar.V_vox_sys{i} = repmat(BodyPar.COG(i,:),size(BodyPar.V_vox_cog{i},1),1)+BodyPar.V_vox_cog{i};
        V_vox_sys = BodyPar.V_vox_sys{i};
        dm = BodyPar.dVol(i)*BodyPar.rho(i);
        Ixx_sys = sum((V_vox_sys(:,2).^2 + V_vox_sys(:,3).^2)*dm);
        Iyy_sys = sum((V_vox_sys(:,1).^2 + V_vox_sys(:,3).^2)*dm);
        Izz_sys = sum((V_vox_sys(:,1).^2 + V_vox_sys(:,2).^2)*dm);
        Ixy_sys = -sum((V_vox_sys(:,1).*V_vox_sys(:,2))*dm);
        Ixz_sys = -sum((V_vox_sys(:,1).*V_vox_sys(:,3))*dm);
        Iyz_sys = -sum((V_vox_sys(:,2).*V_vox_sys(:,3))*dm);
        BodyPar.Isys = BodyPar.Isys + [Ixx_sys Ixy_sys Ixz_sys
                                       Ixy_sys Iyy_sys Iyz_sys
                                       Ixz_sys Iyz_sys Izz_sys];
    end
    if ~exist('opt.cordLengthRef') || isempty(opt.cordLengthRef)
        opt.cordLengthRef = max(V(:,1))-min(V(:,1));
    end
    
end


if isfield(opt,'autoSREFFlag') && opt.autoSREFFlag == 0
    opt.SREF = opt.SREF;
    opt.LREF = opt.LREF;
end

if isfield(opt,'MESH') && ~isempty(opt.MESH)
    Fstart = opt.MESH{1};
    Vstart = opt.MESH{2};
    F = opt.MESH{1};
    V = opt.MESH{2};
    if opt.ThermoFlag == 1
        rNverts = opt.MESH{4};
        rNfaces = opt.MESH{5};
        N_verts = opt.MESH{6};

        if opt.ablationFlag == 1
            Tw = opt.MESH{7};
        end
    end
    Afaces = Area(F,V);
    if ~isfield(opt,'cordLengthRef') || isempty(opt.cordLengthRef)
        opt.cordLengthRef = max(V(:,1))-min(V(:,1));
    end
end

if opt.voxFlag == 1 && isfield(opt,'MESH') && ~isempty(opt.MESH)
    % computing the object properties from the system reference frame
    for i = 1:Nbodies
        [BodyPar.V_vox_cog{i}, BodyPar.COG(i,:), BodyPar.I{i}, BodyPar.dVol(i,1), BodyPar.mass(i,1)] = ...
                        Run_Voxelise(BodyPar, BodyPar.F{i}, BodyPar.V_sys{i}, opt.NnormGrid, BodyPar.rho(i), BodyPar.ShellFlag(i), BodyPar.shellThickness(i), i);
    end
    BodyPar.V_cog{i} = BodyPar.V_sys{i} - repmat(BodyPar.COG(i,:), length(BodyPar.V_sys{i}),1);
    BodyPar.InMass = sum(BodyPar.mass);
    BodyPar.Isys = zeros(3,3);
    BodyPar.COGsys = sum(BodyPar.COG.*repmat(BodyPar.mass,1,3),1)./sum(BodyPar.mass);
    BodyPar.COGsys_hist(size(BodyPar.COGsys_hist,1)+1,:)  = BodyPar.COGsys;
    BodyPar.COG = BodyPar.COG-repmat(BodyPar.COGsys, Nbodies,1);
    BodyPar.COG_hist{length(BodyPar.COG_hist)+1,1}  = BodyPar.COG;
    for i = 1:Nbodies
        BodyPar.V_sys{i} = BodyPar.V_sys{i}-repmat(BodyPar.COGsys, length(BodyPar.V_sys{i}),1);
        BodyPar.V_vox_sys{i} = repmat(BodyPar.COG(i,:),size(BodyPar.V_vox_cog{i},1),1)+BodyPar.V_vox_cog{i};
        V_vox_sys = BodyPar.V_vox_sys{i};
        dm = BodyPar.dVol(i)*BodyPar.rho(i);
        Ixx_sys = sum((V_vox_sys(:,2).^2 + V_vox_sys(:,3).^2)*dm);
        Iyy_sys = sum((V_vox_sys(:,1).^2 + V_vox_sys(:,3).^2)*dm);
        Izz_sys = sum((V_vox_sys(:,1).^2 + V_vox_sys(:,2).^2)*dm);
        Ixy_sys = -sum((V_vox_sys(:,1).*V_vox_sys(:,2))*dm);
        Ixz_sys = -sum((V_vox_sys(:,1).*V_vox_sys(:,3))*dm);
        Iyz_sys = -sum((V_vox_sys(:,2).*V_vox_sys(:,3))*dm);
        BodyPar.Isys = BodyPar.Isys + [Ixx_sys Ixy_sys Ixz_sys
                                       Ixy_sys Iyy_sys Iyz_sys
                                       Ixz_sys Iyz_sys Izz_sys];
    end
    if ~isfield(opt,'cordLengthRef') || isempty(opt.cordLengthRef)
        opt.cordLengthRef = max(V(:,1))-min(V(:,1));
    end
end

if opt.voxFlag == 1 || isempty(opt.MESH)
    V = V-repmat(BodyPar.COGsys,length(V),1);    
end
CofR_offset = [min(V(:,1)) 0 0];
CofR = [0 0 0]+CofR_offset; %BodyPar.COGsys;

faceBaricenter(:,1) = sum(reshape(V(F(:,1:3),1),[],3),2)/3;
faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3),2)/3;
faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3),2)/3;


for i = 1:Nbodies
    % Computing the local object barycenter
    BodiesCOG_matr_F{i} = repmat(BodyPar.COG(i,:),BodyPar.LimtsF(i,2)-BodyPar.LimtsF(i,1)+1,1);

    if ~isempty(opt.fixed_COG)
        warning('change the fixed COG ');
        BodiesCOG_matr_F{i} = opt.fixed_COG;
    end
    
%     % Normals Computation Alternative    
%     opt.MESH.faces = F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:);
%     opt.MESH.vertices = V;
%     
%     % Corrects also the winding direction
%     [N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)] = COMPUTE_opt.MESH_normals(opt.MESH);

%     % Normals Computation
    e0=V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),3),:)-V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),2),:);
    e1=V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),1),:)-V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),3),:);
    Ntemp=cross(e1,e0);
    N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)=normr(Ntemp); 
    
    % checking if the normals have been computed correctly, flat face normals may have an inverse sign
    faceNrm = normr(faceBaricenter(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)-BodiesCOG_matr_F{i});
    normalsCheck = dot(faceNrm,N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:),2);
   
%     opt.MESH.faces = F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:);
%     opt.MESH.vertices = V;
%     [N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)] = COMPUTE_opt.MESH_normals(opt.MESH);    
    
    % %%%%%%%% TO BE REACTIVATED FOR OTHER GEOMETRIES %%%%%%%%%%%
    %         Changing the winding direction
    if sum(normalsCheck<0) > length(N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:))/2
        Ftemp = F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),3);
        F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),3) = F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),2);
        F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),2) = Ftemp;
        
        e0=V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),3),:)-V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),2),:);
        e1=V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),1),:)-V(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),3),:);
        Ntemp=cross(e1,e0);
        N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)=normr(Ntemp);         
        
%         opt.MESH.faces = F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:);
%         opt.MESH.vertices = V;
%         [N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)] = COMPUTE_opt.MESH_normals(opt.MESH);
    end
    
    % %%%%%%%% opt.MESH coherence check for simple geometries %%%%%%%%%%%
    %         Changing the winding direction
    if BodyPar.SimpleGeom(i) == 1
        if sign(sum(normalsCheck)) == 1
            vertToChange = find(normalsCheck < 0);
        else
            vertToChange = find(normalsCheck > 0);
        end
        Ftemp = F(BodyPar.LimtsF(i,1)-1+vertToChange,3);
        F(BodyPar.LimtsF(i,1)-1+vertToChange,3) = F(BodyPar.LimtsF(i,1)-1+vertToChange,2);
        F(BodyPar.LimtsF(i,1)-1+vertToChange,2) = Ftemp;
        
        % Recomputing the normals
        e0=V(F(BodyPar.LimtsF(i,1)-1+vertToChange,3),:)-V(F(BodyPar.LimtsF(i,1)-1+vertToChange,2),:);
        e1=V(F(BodyPar.LimtsF(i,1)-1+vertToChange,1),:)-V(F(BodyPar.LimtsF(i,1)-1+vertToChange,3),:);
        Ntemp=cross(e1,e0);
        N(BodyPar.LimtsF(i,1)-1+vertToChange,:)=normr(Ntemp);                  
    end    
    
    BodyPar.N{i} = N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:);
end


if opt.checks == 1
    fig = figure();
    fig.Position = [100, 25, 950, 700];
    trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); axis('equal'); hold on;
    plot3(BodyPar.COG(:,1),BodyPar.COG(:,2),BodyPar.COG(:,3),'or','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','r')
    plot3(0,0,0,'sb','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','b')
    quiver3(faceBaricenter(:,1),faceBaricenter(:,2),faceBaricenter(:,3), N(:,1), N(:,2), N(:,3),8,'b');
    axis([min(V(:,1)) max(V(:,1)) (min(V(:,2))+max(V(:,2)))/2 max(V(:,2)) min(V(:,3)) max(V(:,3))]*1.4)
%     axis([min(V(:,1)) max(V(:,1)) min(V(:,2)) (min(V(:,2))+max(V(:,2)))/2  min(V(:,3)) max(V(:,3))]*1.2)
    disp('Normals Check');
    keyboard
    
% %     % If some faces have a wrong winding use the following to change them (for each body):
% %     indToChange = find(normalsCheck<0);
% %         Ftemp = F(indToChange,3);
% %         F(indToChange,3) = F(indToChange,2);
% %         F(indToChange,2) = Ftemp;
% %         
% %         opt.MESH.faces = F;
% %         opt.MESH.vertices = V;
% %         N = COMPUTE_opt.MESH_normals(opt.MESH);            
end

%%

if isempty(opt.MESH)
    opt.MESH{1} = F;
    opt.MESH{2} = V; 
    opt.MESH{3} = N;    
    Fstart = F;
    Vstart = V;
end

if (opt.autoSREFFlag == 1 && ((isempty(opt.SREF) || opt.SREF == 0 ) || (isempty(opt.LREF) || opt.LREF == 0 ))) ...
        || ((isempty(opt.SREF) || opt.SREF == 0 ) || (isempty(opt.LREF) || opt.LREF == 0 ))
    [Fopt.SREF, ~, Nopt.SREF] = BackFaceCulling_v2_2(0, 0, 0, F, V, N, 1500, 0 , opt.bfcFlag, 0);    
    Nf = cross(repmat([1 0 0],length(Nopt.SREF),1),Nopt.SREF());
    Theta = asind(sqrt(Nf(:,1).^2+Nf(:,2).^2+Nf(:,3).^2)); 
    Afront = Area(Fopt.SREF,V);
    if isempty(opt.SREF) || opt.SREF == 0 
        opt.SREF = sum(Afront.*sind(90-Theta));  
        opt.SREF = opt.SREF;
    else
        opt.SREF = opt.SREF;
    end
    if isempty(opt.LREF) || opt.LREF == 0 
        opt.LREF = max([max(V(:,2))-min(V(:,2)) max(V(:,3))-min(V(:,3))]);  
        opt.LREF = opt.LREF;
    else
        opt.LREF = opt.LREF;
    end
else
    opt.SREF = opt.SREF;
    opt.LREF = opt.LREF;
end

% if not provided, use a reference maximum flat radius
if ~isfield(opt,'RmaxRef') || isempty(opt.RmaxRef)
    opt.RmaxRef = sqrt(opt.SREF/pi)*2;
end
Aratio = opt.SREF/sum(Afaces);

%% Computing the local radius
if ~exist('rNfaces','var')
    if opt.ThermoFlag ==1 
        if opt.NsmoothRadius ==0
            SmoothFlag = 0;
        else
            SmoothFlag = 1;
        end
    for i = 1:Nbodies
        if i == 1
            [FVNN, RadiiV_Vs_F_Fs] = STLcurvature(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:), ...
                V(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:), N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:),...
                opt.RmaxRef/M2M_RR, opt.RmaxRef, {SmoothFlag, round(opt.NsmoothRadius*length(F)/1000), 'e', opt.FlatEdge, 2}, 0);
        else
            [FVNN, RadiiV_Vs_F_Fs] = STLcurvature(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)-BodyPar.LimtsV(i-1,2),...
                V(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:), N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:),...
                opt.RmaxRef/M2M_RR, opt.RmaxRef, {SmoothFlag, round(opt.NsmoothRadius*length(F)/1000), 'e', opt.FlatEdge, 2}, 0);
        end
        N_verts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = FVNN{4};
        rNverts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = RadiiV_Vs_F_Fs{2};
        rNfaces(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:) = RadiiV_Vs_F_Fs{4};
        opt.MESH{4}(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = rNverts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:);
        opt.MESH{5}(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:) = rNfaces(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:);
        opt.MESH{6}(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = N_verts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:);
    end
    end
elseif exist('rNfaces','var') && opt.radiusUpdateFlag == 1
    if opt.NsmoothRadius ==0
        SmoothFlag = 0;
    else
        SmoothFlag = 1;
    end
    for i = 1:Nbodies
        if i == 1
            [FVNN, RadiiV_Vs_F_Fs] = STLcurvature(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:), ...
                V(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:), N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:),...
                opt.RmaxRef/M2M_RR, opt.RmaxRef, {SmoothFlag, round(opt.NsmoothRadius*length(F)/1000), 'e', opt.FlatEdge, 2}, 0);
        else
            [FVNN, RadiiV_Vs_F_Fs] = STLcurvature(F(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:)-BodyPar.LimtsV(i-1,2),...
                V(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:), N(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:), opt.RmaxRef/M2M_RR,...
                opt.RmaxRef, {SmoothFlag, round(opt.NsmoothRadius*length(F)/1000), 'e', opt.FlatEdge, 2}, 0);
        end
        N_verts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = FVNN{4};
        rNverts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = RadiiV_Vs_F_Fs{2};
        rNfaces(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:) = RadiiV_Vs_F_Fs{4};
        opt.MESH{4}(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = rNverts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:);
        opt.MESH{5}(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:) = rNfaces(BodyPar.LimtsF(i,1):BodyPar.LimtsF(i,2),:);
        opt.MESH{6}(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:) = N_verts(BodyPar.LimtsV(i,1):BodyPar.LimtsV(i,2),:);
    end
end

    if opt.checks == 1 && opt.ThermoFlag == 1
        fig = figure();
        fig.Position = [100, 50, 950, 700];
        trisurf(F,V(:,1),V(:,2),V(:,3),rNfaces,'edgealpha',0.15); axis('equal')
        colorbar; set(gca, 'CLim', [min(rNfaces) max(rNfaces)]);
        disp('Check the local radius quality')
        keyboard
    end

%% Facets visibility detection

er = 1;
AreaError = 1;
areaBfcOld = 0;
Nit=0;
% Area based BFC convergence and refinement loop
while AreaError(er) >= AreaConvError && n_pix_temp(er) < maxPixRes+1
    tsimbfc = tic;
    if length(BodyPar.ablFaceID_sys_global{end}) > 15
        Ftemp  = F(BodyPar.survingFaceID_sys_global{end},:);
        Ntemp = N(BodyPar.survingFaceID_sys_global{end},:);
        [bfcFVN{1}, ~, bfcFVN{3}] = BackFaceCulling_v2_2(AofA, SS, Roll, Ftemp, V, Ntemp, n_pix_temp(er), 0, opt.bfcFlag, 1);
    else    
        [bfcFVN{1}, ~, bfcFVN{3}] = BackFaceCulling_v2_2(AofA, SS, Roll, F, V, N, n_pix_temp(er), 0, opt.bfcFlag, 0);
    end
    tsimbfc = toc(tsimbfc);
    %             nFaceOld(er+1) = length(bfcFVN{1});
    %             nFaceError(er+1) = abs((nFaceOld(er+1)-nFaceOld(er))./nFaceOld(er+1));
    areaBfcOld(er+1) = sum(Area(bfcFVN{1},V));
    AreaError(er+1) = abs((areaBfcOld(er+1)-areaBfcOld(er))./areaBfcOld(er+1));
%     disp(['Total Area detected by BFC: ',num2str(areaBfcOld(er+1))])
%     disp(['Total Area convergence error with ',num2str(n_pix_temp(er)),' pixels: ',num2str(AreaError(er+1))])
    n_pix_temp(er+1) = n_pix_temp(er)+opt.nPixStep;
    er = er+1;
    Nit = Nit+1;
    
    if AreaError(er) <= AreaConvError || n_pix_temp(er) >= maxPixRes+1
        shadowRefFlag = 0;
    end

    if shadowRefFlag == 1
        
%%         % Not Tested in this version
%         
%         %finding the shadowed-bfc border indexes
%         [~,IB,preBFCind] = intersect(bfcFVN{1},F,'rows','stable');
%         Fbfc = F(preBFCind,:);
%         Nbfc = bfcFVN{3}(IB,:);
%         shadowInd = setdiff([1:length(F)]',preBFCind);
%         Fshadow = F(shadowInd,:);
%         Nfaces = length(F);
%         
%         [refIndex,~,~] = intersect(reshape(Fbfc,[],1),reshape(Fshadow,[],1));
%         [refIndex] = find(ismember(F,refIndex));
%         refIndex(refIndex > Nfaces & refIndex <= Nfaces*2) = refIndex(refIndex > Nfaces & refIndex <= Nfaces*2)-Nfaces;
%         refIndex(refIndex > Nfaces*2) = refIndex(refIndex > Nfaces*2)-Nfaces*2;
%         refIndex = unique(refIndex,'stable');
%         
%         % Double bordering -- - shock impingement purposes
%         % % %                 [refIndex,~,~] = intersect(reshape(F,[],1),reshape(F(refIndex,:),[],1));
%         % % %                 [refIndex] = find(ismember(F,refIndex));
%         % % %                 refIndex(refIndex > Nfaces & refIndex <= Nfaces*2) = refIndex(refIndex > Nfaces & refIndex <= Nfaces*2)-Nfaces;
%         % % %                 refIndex(refIndex > Nfaces*2) = refIndex(refIndex > Nfaces*2)-Nfaces*2;
%         
%                     if opt.checks == 1
%                         fig = figure();
%                         fig.Position = [100, 100, 950, 800];
%                         trisurf(Fbfc,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
%                         trisurf(F(shadowInd,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]);
%                         trisurf(F(refIndex,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0 0 1]);
%                         trisurf(F(shadowInd,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]);
%                         disp('Shadowing Refinement check Triangular opt.MESH Check');
%                         keyboard
%                     end
%         [F, V] = cleanVerts(F, V);
%         [F, V, N] = refineSTL_tri_v2_2(F, V, 'both',[Amax Lmax],length(F)+length(refIndex)*Nref,10,{'on',refIndex});
%         
% %     fig = figure();
% %     fig.Position = [100, 100, 950, 800];
% %     trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); axis('equal'); hold on;
% %     faceBaricenter = sum(reshape(V(F(:,1:3),1),[],3),2)/3;
% %     faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3),2)/3;
% %     faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3),2)/3;
% %     fq1 = quiver3(faceBaricenter(:,1),faceBaricenter(:,2),faceBaricenter(:,3), N(:,1), N(:,2), N(:,3),3,'r');        
%%
    end
end

bfcFVN{2} = V;


Afaces = Area(F,V);
faceBaricenter = sum(reshape(V(F(:,1:3),1),[],3),2)/3;
faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3),2)/3;
faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3),2)/3;

[~,IB,preBFCind] = intersect(bfcFVN{1},F,'rows','stable');
Fbfc = F(preBFCind,:);
Nbfc = bfcFVN{3}(IB,:);
if length(Tw) > 1
    Tw = Tw(preBFCind);
end
Nfaces = length(Nbfc);
faceBarBFC = faceBaricenter(preBFCind,:);
Abfc = Afaces(preBFCind);
Atot = sum(Abfc);

% figure()
% trisurf(F,V(:,1),V(:,2),V(:,3),Afaces,'edgealpha',0.15);
% trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Abfc,'edgealpha',0.15);

if opt.checks == 1
    shadowInd = setdiff([1:length(F)]',preBFCind);
    fig = figure();
    fig.Position = [100, 100, 950, 800];
    trisurf(bfcFVN{1},V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
    trisurf(F(shadowInd,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]);
    Labels = {'x','y','z'};
    Title = {[STLnames{1}(1:end-4),' Shadowing Refinement']};
    Annotation = {['$\alpha = ',num2str(AofA),'deg$ and $\beta = ',num2str(SS),'deg$']};
    Legend = {'Non-Shadowed','Shadowed'}; %,'FOSTRAD $A_{ratio} = 0.10$','DSMC $A_{ratio} = 0.10$', ...
    stdPlot(Labels, Title, Legend, Annotation)
    disp('Normals Check');
    keyboard
end


%% Computing flow incidence angle and applying aerodynamic corrective factors
bb=1;
ROT = Rotation_v1(SS,AofA,Roll);
B2WA = ROT';
%Body to Wind rotation matrix
Vinfi = (ROT * [Vinf(1); 0; 0])';
Vinfni = Vinfi./Vinf(1);

% computing the correction factor
if opt.CF_Flag == 0
    dragCF = 1;
elseif opt.CF_Flag == 1
    dragCF = FostradToDSMC_CD_parallelepiped(SS,AofA,Aratio);
elseif opt.CF_Flag == 2
    eqAngle = asind(sqrt(Vinfni(2)^2+Vinfni(3)^2));
    dragCF = FostradToDSMC_CD_Cylinder(eqAngle, Aratio);
end

Nf = cross(repmat(Vinfni,Nfaces,1),Nbfc);
Nf = sqrt(Nf(:,1).^2+Nf(:,2).^2+Nf(:,3).^2);
Nf(Nf>1) = 1;
Theta = asind(Nf); 
    
if opt.ThermoFlag ==1
    minTheta = min(Theta);
    % finding the stagnation point/points
    opt.stagFaceIndvect = find(Theta == minTheta);

    Aview = Abfc.*sind(90-Theta);
    AtotView = sum(Aview);
    rNmax = 2*sqrt(AtotView/pi());
    rNverts(rNverts > rNmax) = rNmax;

% %     % necessary if the shadowing refinement has been used
% %         if shadowRefFlag == 1
% %             tremap = tic;
% %             rNvertsNew = MESHRemapSmooth(Fstart, Vstart, V, rNverts, 1);
% %             tremap = toc(tremap);
% % %             disp(['Remapping local radius on the refined opt.MESH: ',num2str(tremap),'s'])
% %             rNfaces = verts2face(F, V, rNvertsNew);
% %         end         
    rN = rNfaces(preBFCind); 
    
    if opt.checks == 1
        fig = figure();
        fig.Position = [100, 50, 950, 700];
        trisurf(F,V(:,1),V(:,2),V(:,3),rNfaces,'edgealpha',0.15); hold on; axis('equal')
        trisurf(Fbfc,V(:,1),V(:,2),V(:,3),rN,'edgealpha',0.15);
        colorbar; set(gca, 'CLim', [min(rNfaces) max(rNfaces)]);
        disp('Radius Quality Check');
        keyboard
    end
end


%% Computing Atmospheric Properties

if size(altitude,2) > size(altitude,1)
    altitude = altitude';
end

% Rearranging the database inputs to match FOSTRAD inputs
%He O N2 O2 Ar H N [1/m^3] gas densities
rhoi = [NRLMSISE00(Hcont+1:Hfmf+1,6) NRLMSISE00(Hcont+1:Hfmf+1,2:4) NRLMSISE00(Hcont+1:Hfmf+1,7:9)]*1E6;
Ti = NRLMSISE00(Hcont+1:Hfmf+1,5);

%inputs corresponding to the input altitudes interpolated with the
%continuum and FMF reference
hDB = NRLMSISE00(Hcont+1:Hfmf+1,1); %altitudes used for the sample points interpolation on the database
h= unique([linspace(min(hDB),max(hDB),1e4)'; altitude]);
if length(hDB) > 1
    Ti = interp1(hDB,Ti(:,1),h,'pchip');
    for i=1:7
        rhotemp(:,i) = interp1(hDB,rhoi(:,i),h,'pchip');
    end
    rhoi = rhotemp;
end


% % % opt.Env = MixConstants([rhoi(:,2:5), rhoi(:,7), rhoi(:,1), rhoi(:,6)], Ti(:), VHSflag);
if ~isfield(opt,'Env') || isempty(opt.Env) || size(rhoi,1) ~= size(opt.Env,1)
    opt.Env = MixConstants_v([rhoi(:,2:5), rhoi(:,7), rhoi(:,1), rhoi(:,6)], Ti(:), VHSflag);
end
%[C0' C1' C2' C3' RHO' P' Rmix' Cpmix' GAMMA' MMean' OMEGA' Vis'];
C = [opt.Env(:,1:4)]; % constants for adimensionals
rho = opt.Env(:,5); %kg/m3
P = opt.Env(:,6); %Pa
R = opt.Env(:,7); %j/molK
cp = opt.Env(:,8); %j/kgK
gamma = opt.Env(:,9); 
m = opt.Env(:,10); %kg/mol
omega  = opt.Env(:,11); 
muEC  = opt.Env(:,12); %pa.s Enskog-Chapman transitional/FMF
muSu = opt.Env(:,13); %pa.s Sutherlands Continuum
mu = zeros(length(h),1);  %initialising
% Vcir = sqrt(398600./(h+6378))*1000;

% Speed Ratio or Velocity Check or orbital velocity
if Vinf == 0 % Use Circular Orbit velocity Vinf = F(H)
    disp('Using Circular orbit velocity')
    Vinf = sqrt(398600./(h+6378))*1000;
    SRflag = 0;
elseif Vinf < 50 && Vinf > 0 %use constant speed ratio
%     disp('Using Constant speed ratio')
    SR = ones(length(h),1) .* Vinf;
    Vinf = Vinf .* C(:,4) .* sqrt(Ti);  %applying the Speed ratio formulation
    SRflag = 1;
elseif Vinf >= 50       % use the provided velocity
%     disp('Using provided constant velocity')
    Vinf = ones(length(h),1) .* Vinf;
    SRflag = 0;
end

%Adimensionals = [ Re, SR, Mach, Kn, Tratio]
adims = Adimensionals_v1(opt.Env(:,1:5), Ti(:), Tw, Vinf(:), opt.LREF);
if exist('SR','var') == 0
    SR = adims(:,2);
end
Ma = adims(:,3);
Kn = adims(:,4);

% Building the viscosities in the different regimes
MuBoundMin = find(Kn >= 0.08, 1,'first');
MuBoundMax = find(Kn >= 8, 1,'first');
bridge_Ind = [MuBoundMin-3:MuBoundMin, MuBoundMax:MuBoundMax+3];
mu_bridge = pchip(Kn(bridge_Ind), [muSu(bridge_Ind(1:4))' muEC(bridge_Ind(5:8))'], Kn(MuBoundMin:MuBoundMax));

mu(1:MuBoundMin) = muSu(1:MuBoundMin);
mu(MuBoundMin:MuBoundMax) = mu_bridge;
mu(MuBoundMax:end) = muEC(MuBoundMax:end);

%   Thermal Conductivity Coefficient
k = 2.64638e-3*Ti.^1.5./(Ti+245*10.^(-12./Ti));
Pr = mu.*cp./k;


% Free Stream Stagnation Pressure and Temperature
% immediately after the shock wave
% P01 = P .* (2 * gamma .* Ma.^2 - (gamma-1)) ./ (gamma+1);

% at the stagnation point
P02 = P .*(0.5*(gamma+1).*Ma.^2).^(gamma./(gamma-1)) ... 
             .* ((gamma+1)./(2*gamma.*Ma.^2 - (gamma-1))).^(1./(gamma-1));

% temperature at the stagnation point  
T0s = Ti .* (1 + 0.5 .* (gamma-1) .* Ma.^2);

mu_T0s = ViscosityWall([rhoi(:,2:5), rhoi(:,7), rhoi(:,1), rhoi(:,6)], T0s.*ones(length(rhoi(:,2:5)),1));   % Stcagnation point viscosity through Sutherland's law (continuum) 
% mu_T0 = mu.* (T0./Ti).^omega;  % viscosity via the Chapman-Enskog eq. not valid in the continuum
% rhow = P02./R./Tw;  % Density at the wall boundary layer

% hw = cp .* Tw;
rhos = P02./R./T0s;  % Density at the stagnation point
h0s = cp .* T0s;
Re0norm = rho .* Vinf .* 1 ./ (mu.*(T0s./Ti).^omega);   % Reynolds for SCARAB formulation "SCARAB - A Multi-disciplinary code for destruction analysis..."

% Stagnation point pressure computation
Cpmax=(2./(gamma.*Ma.^2)).*((P02./P-1)); %Cp maximum


%% 
% subdividing the interval according to knudsen; defininig points near and far from the transition points
KnDiv = [limKn_inf*0.99 limKn_inf_heat*0.99 limKn_sup*1.01];

KnRef = repmat(Kn,1,length(KnDiv));
KnDiv = repmat(KnDiv,length(Kn),1);
[~, hInd] =  min(abs(KnRef-KnDiv));

hIndReq = find(sum(repmat(h,1,length(altitude)) == repmat(altitude',length(h),1),2) == 1);
contFlag = 0;
FMFflag = 0;
if  all(Kn(hIndReq) >= limKn_sup) 
    bridgeFlag = 0;
    FMFflag = 1;
    hInd = hIndReq;
elseif  all(Kn(hIndReq) <= limKn_inf)
    bridgeFlag = 0;
    contFlag = 1;
    hInd = hIndReq;
else
    bridgeFlag = 1; 
    hInd = unique([hInd'; hIndReq]);
end

altitude = h(hInd);
NKnc = length(find(Kn(hInd) <= limKn_inf));
NKntrans = length(find(Kn(hInd) > limKn_inf & Kn(hInd) < limKn_sup ));

NKnFMF = length(find(Kn(hInd) >= limKn_sup ));

Knc_h_ind = find(Kn(hInd) <= limKn_inf_heat);
NKnc_h = length(Knc_h_ind);
Knc_h = Kn(hInd(Knc_h_ind));

Kntrans_h_ind = find(Kn(hInd) > limKn_inf_heat & Kn(hInd) < limKn_sup);
NKntrans_h = length(Kntrans_h_ind);
Kn_trans_h = Kn(hInd(Kntrans_h_ind));

%%  Parameter Initialization

Cpc = zeros(Nfaces,1);

if opt.ThermoFlag == 1
    Stc = zeros(Nfaces,1);
    Stc_FMF = zeros(Nfaces, NKnFMF);
end

Nalt = length(hInd);
if opt.ThermoFlag == 1
    
    % Used for performing the constant radius tests
    if opt.RN_const_flag == 1 && ~isempty(opt.rN_ref)
        rN(:) = opt.rN_ref;
    elseif opt.RN_const_flag == 1
        rN(:) = opt.RmaxRef;
    end
    mu_w = ViscosityWall_v1([rhoi(hInd,2:5), rhoi(hInd,7), rhoi(hInd,1), rhoi(hInd,6)], Tw);   %Function for computing viscosity using sutherlands formulation for the boundary layer
    dudx = repmat((1./rN),1,Nalt) .* repmat((sqrt(2.*(P02(hInd)-P(hInd))./rhos(hInd)))',Nfaces,1);
    Re0 = repmat(2*rN,1,Nalt) .* repmat(Re0norm(hInd)',Nfaces,1);   % Defined with the Diameter for SCARAB formulation

    % stanton number constant: Stc = Q./StConst
        
% %         % This stanton number formulation is valid only for the continuum regime and not for the FMF
% %         StConst = repmat(rho(hInd)' .* Vinf(hInd)' .* cp(hInd)', Nfaces, 1) .* (repmat(T0s(hInd)', Nfaces, 1) - repmat(Tw,1,Nalt)); % ~= CH_const

        % The Heat transfer coefficient is more general and the bridging function have been calibrated on this parameter; although their values should be very similar
        CH_const = rho(hInd).*Vinf(hInd).^3 / 2;
        StConst = repmat(CH_const',Nfaces,1);

    StConst(StConst <= 0.05) = 0.05;  % neglecting the cooling of the surface
end

% CH_const = 2./(rho.*Vinf.^3);

%He O N2 O2 Ar H N [1/m^3] gas densities

%% Computation at different altitudes
    a=1;
    bbb=1;
    
    Delta = 90 - Theta;   % flow inclination
    r = faceBarBFC - repmat(CofR,Nfaces,1);    %Moment Arm
    
    Vinfni = repmat(Vinfni,Nfaces,1);

    % If only one altitude has been requested then the other altitude are used only for building the 
    % bridging functions, therefore the speed ratio is assumed to be equal to the one of the requested
    % altitude throughout all regimes, in this way the briding can be computed correctly
    if length(hIndReq) == 1
        SRi = SR(hIndReq)*ones(Nalt,1);
        Vinfi = Vinf(hIndReq)*ones(Nalt,1);
    else
        SRi= SR(hInd);
        Vinfi = Vinf(hInd);
    end   

for hh=1:Nalt %Km
    H = hInd(hh);

    % Parameter Re-Initialization
    Cc = zeros(1,3);
    Cfm = zeros(1,3);
    Mc = zeros(1,3);
    Mfm = zeros(1,3);
    rhowi = P02(H)/R(H)./Tw;   
    
    %% ******************* AERODYNAMICS MODULE ********************** %%
    if opt.AeroFlag ==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% CONTINUUM AERODYNAMICS %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  Kn(H) <= limKn_inf
            Cpc = Cpmax(H) * (sind(Delta)).^2;
            Cpc(Delta <= 0,1) = 0;            
            Cc = repmat(Cpc,1,3) .* -Nbfc .* repmat(Abfc,1,3);
            Mc = cross(r,Cc);        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%%%%%% FREE MOLECULAR AERODYNAMICS %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif Kn(H) >= limKn_sup
            % Free Molecular Heat Transfer for an inclined flat plate - Ref: Hypersonic Flow Theory - Wallace Hayes, page 402
            t = (Nbfc .* repmat(dot(Vinfni,Nbfc,2),1,3) - Vinfni) ./ repmat(sqrt(1 - dot(Vinfni,Nbfc,2).^2),1,3);
            Cpfm1 = exp(-(SRi(hh)*sind(Delta)).^2) .* ((2-SN)*SRi(hh)*sind(Delta)./sqrt(pi) +...
                SN.*sqrt(Tw/Ti(H))/2);
            Cpfm2 = (1+erf(SRi(hh)*sind(Delta))) .* ((2-SN) .* ((SRi(hh)*sind(Delta)).^2 + 0.5) +...
                0.5*SN*SRi(hh)*sind(Delta).*sqrt(pi*Tw/Ti(H)));
            Cpfm = (Cpfm1+Cpfm2)/SRi(hh)^2;
            Ctfm = -(ST*cosd(Delta)/SRi(hh)/sqrt(pi)) .* (exp(-(SRi(hh)*sind(Delta)).^2) + sqrt(pi) *...
                SRi(hh) *sind(Delta).*(1+erf(SRi(hh)*sind(Delta))));
            t(isnan(t)) = 0;
            
            Cpnfm = repmat(Cpfm,1,3) .* -Nbfc .* repmat(Abfc,1,3);
            Cttfm = repmat(Ctfm,1,3) .* t .* repmat(Abfc,1,3);
            Cfm =  Cpnfm + Cttfm;
            Mfm =  cross(r,Cfm);
        end
    end %end of Aerodynamics module
    
    %% ********************** Aero-ThermoDynamics Module *********************%%
    if opt.ThermoFlag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% CONTINUUM AERODYNAMICS HEATING %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Kn(H) < limKn_inf_heat
            switch lower(opt.AeroThModel)
                case {'sc'}
                    Stc(:,hh) = 2.1./sqrt(Re0(:,hh));             %Stagnation Point Stcanton Number, SCARAB Formulation
                    Stc(:,hh) = Stc(:,hh) .* (0.1+0.9*cosd(Theta)); % Stanton with SCARAB heat formulation and Lees distribution
                case {'fr'} % fully catalytic
                    %               Fay-Riddel "theory of stagnation point heat transfer in dissociated air"
                    Stc(:,hh) = 0.763.*(Pr(H).^-0.6).*((rhos(H).*mu_T0s(H)).^0.4).*((rhowi.*mu_w(hh)).^0.1).*sqrt(dudx(:,hh)).*(h0s(H)-cp(H).*Tw)./StConst(:,hh);
                    Stc(:,hh) = Stc(:,hh) .* (0.1+0.9*cosd(Theta)); % Stanton with Fay-Riddel heat formulation and Lees distribution
                case {'krd'} % super catalytic
                    %               R. W. Detra, N. H. Kemp, F. R. Riddell, Addendum to heat transfer to satellite vehicles reentering the atmosphere, 
                                    % Jet Propulsion, vol. 27, (1957)   -   highly dependant on the ambient density
                    Stc(:,hh) = 110.34*1e6./sqrt(rN(:)).*sqrt(rho(H)./rhosl).*(Vinfi(hh)./7925).^3.15.*(h0s(H)-cp(H).*Tw)./(h0s(H)-cp(H)*300)./StConst(:,hh);  
                    Stc(:,hh) = Stc(:,hh) .* (0.1+0.9*cosd(Theta));
                case {'vd'}
                    % Van Driest - non catalytic
                    Stc(:,hh) = 0.763*(Pr(H)^-0.6)*(rhos(H)*mu_T0s(H))^0.5.*sqrt(dudx(:,hh)).*(h0s(H)-cp(H).*Tw)./StConst(:,hh);
                    Stc(:,hh) = Stc(:,hh) .* (0.1+0.9*cosd(Theta));
            end
%         Stc = Stc(i,hh) .* (cosd(Theta/2)^5.27); % Stanton using Kemp Rose Detra distribution (as in ORSAT)
        Stc(Stc(:,hh)<0, hh) = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% FMF THERMODYNAMICS HEATING %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif Kn(H) >= limKn_sup
            Qstfm = 0.5 * opt.AccCoeff .* rho(H) .* Vinfi(hh).^3; % acc coefficient = 0.93;
            % Free Molecular Heat Transfer for an inclined flat plate - Ref: Hypersonic Flow Theory - Wallace Hayes, page 403
            Q_fm = (Qstfm/(SRi(hh)^3)/2/sqrt(pi))*((repmat(SRi(hh)^2+gamma(H)/(gamma(H)-1),Nfaces,1)-((gamma(H)+1).*Tw/2/(gamma(H)-1)/Ti(H))).*...
                (exp(-(SRi(hh)*sind(Delta)).^2)+sqrt(pi)*SRi(hh)*sind(Delta).*(1+erf(SRi(hh)*sind(Delta))))-0.5.*exp(-(SRi(hh)*sind(Delta)).^2));
                        
            Stc_FMF(:,hh-NKnc_h-NKntrans_h)= Q_fm(:,1)./StConst(:,hh); %heat transfer coefficient fm            
        end
    end %End of ThermoDynamics Module
    
    bb=bb+1; % Altitude indexing for the transitional heating
    if Kn(H)>=limKn_sup
        Knfmf(a,1)=Kn(H);
        a=a+1;
    elseif Kn(H)<=limKn_inf
        Knc(bbb,1)=Kn(H);
        bbb=bbb+1;
    end
    
    %%  Computation of Aerodynamics bridging parameters
    
    if opt.AeroFlag == 1
        CfmTot = sum(Cfm,1);
        MfmTot = sum(Mfm,1);
        CcTot = sum(Cc,1);
        McTot = sum(Mc,1);
        CcTot = CcTot/opt.SREF;
        CfmTot = CfmTot/opt.SREF;
        McTot = McTot/opt.SREF/max(opt.cordLengthRef);
        MfmTot = MfmTot/opt.SREF/max(opt.cordLengthRef);
        CDc = dot(CcTot, B2WA(1,:));
        CSc = dot(CcTot, B2WA(2,:));
        CLc = dot(CcTot, B2WA(3,:));        
        CSfm = dot(CfmTot, B2WA(2,:));
        CLfm = dot(CfmTot, B2WA(3,:));        
        CDfm = dot(CfmTot, B2WA(1,:));
        Mxc = dot(McTot, B2WA(1,:));
        Myc = dot(McTot, B2WA(2,:));
        Mzc = dot(McTot, B2WA(3,:));
        Mxfm = dot(MfmTot, B2WA(1,:));
        Myfm = dot(MfmTot, B2WA(2,:));
        Mzfm = dot(MfmTot, B2WA(3,:));
        if Kn(H) >= limKn_sup
            CDfm_memo(hh-(NKnc+ NKntrans))=[CDfm']*dragCF;
            CLfm_memo(hh-(NKnc+ NKntrans))=[CLfm'];            
            CSfm_memo(hh-(NKnc+ NKntrans))=[CSfm'];
            CMxfm_memo(hh-(NKnc+ NKntrans))=[Mxfm']*dragCF;
            CMyfm_memo(hh-(NKnc+ NKntrans))=[Myfm']*dragCF;
            CMzfm_memo(hh-(NKnc+ NKntrans))=[Mzfm']*dragCF;
        elseif Kn(H) <= limKn_inf
            CDc_memo(hh)=[CDc'];
            CLc_memo(hh) =[CLc'];            
            CSc_memo(hh)=[CSc'];
            CMxc_memo(hh)=[Mxc'];
            CMyc_memo(hh)=[Myc'];
            CMzc_memo(hh)=[Mzc'];
        end
    end
    
    %     figure();
    %     trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
    %     trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Theta,'edgealpha',0.15); colorbar;
    %     figure();
    %     trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
    %     trisurf(Fbfc,V(:,1),V(:,2),V(:,3),rN,'edgealpha',0.15); colorbar; set(gca, 'CLim', [min(rN) max(rN)]);
    %     figure();
    %     trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
    %     trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Stc,'edgealpha',0.15); colorbar; set(gca, 'CLim', [min(min([Stc Stc])) max(max([Stc Stc]))]);
    %     stdPlot({'x','y','z'},{'Scarab'})
    %     figure();
    %     trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
    %     trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Stc,'edgealpha',0.15); colorbar; set(gca, 'CLim', [min(min([Stc Stc])) max(max([Stc Stc]))])
    %     stdPlot({'x','y','z'},{'Fay-Riddel'})
    
end % altitude analyses ended


if bridgeFlag == 1
        % Setting up the bridging functions
        Kn_trans = Kn(hInd);
        Kn_trans = [Knc(end)*1.00001; Kn_trans(Kn_trans > Knc(end) & Kn_trans < Knfmf(1)); Knfmf(1)*0.99999];
        Kn_trans_h = [Knc_h(end)*1.00001; Kn_trans_h; Knfmf(1)*0.99999];
        Kn_trans_R = (log(Kn_trans)-log(Knc(end)))./(log(Knfmf(1))-log(Knc(end)));



%% $%^ Bridging function Correction and calibration
% used to correct the bridging function to match the DSMC transitional regime$%^
        
        % Correction factor to be applied to the bridging function (erf*bridgeCF)
        BridgeCF = Kn_trans_R./((1+erf(Kn_trans_R*4-2))/2); 
        BridgeCF(BridgeCF>1) =  (BridgeCF(BridgeCF>=1)-1)*CF_ratiolow+1;
        BridgeCF(BridgeCF<1) =  (BridgeCF(BridgeCF<=1)-1)*CF_ratiohigh+1;
        
        AeroBridge = (1+erf(Kn_trans_R*4-2))/2.*BridgeCF;
        AeroBridge = (AeroBridge-min(AeroBridge))/(max(AeroBridge)-min(AeroBridge));


    %% ************************ AERODYNAMICS  ****************************%%%%
    %      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.AeroFlag == 1;
        % fixing the sign of the moments
        CD_transERF = [CDc_memo(end)+(CDfm_memo(1)-CDc_memo(end)).*AeroBridge]';
        CL_transERF = [CLc_memo(end)+(CLfm_memo(1)-CLc_memo(end)).*AeroBridge]';
        CS_transERF = [CSc_memo(end)+(CSfm_memo(1)-CSc_memo(end)).*AeroBridge]';
        
        CMx_transERF = [CMxc_memo(end)+(CMxfm_memo(1)-CMxc_memo(end)).*AeroBridge]';
        CMy_transERF = [CMyc_memo(end)+(CMyfm_memo(1)-CMyc_memo(end)).*AeroBridge]';
        CMz_transERF = [CMzc_memo(end)+(CMzfm_memo(1)-CMzc_memo(end)).*AeroBridge]';
        
        [KnBridge, KnInd] = unique([Knc; Kn_trans(2:end-1); Knfmf],'stable');
        
        CD_Bridge = [CDc_memo CD_transERF(2:end-1) CDfm_memo];
        CD_Bridge = CD_Bridge(KnInd);

        CL_Bridge = [CLc_memo CL_transERF(2:end-1) CLfm_memo];
        CL_Bridge = CL_Bridge(KnInd);

        CS_Bridge = [CSc_memo CS_transERF(2:end-1) CSfm_memo];
        CS_Bridge = CS_Bridge(KnInd);        

        CMx_Bridge = [CMxc_memo CMx_transERF(2:end-1) CMxfm_memo];
        CMx_Bridge = CMx_Bridge(KnInd);

        CMy_Bridge = [CMyc_memo CMy_transERF(2:end-1) CMyfm_memo];
        CMy_Bridge = CMy_Bridge(KnInd);
        
        CMz_Bridge = [CMzc_memo CMz_transERF(2:end-1) CMzfm_memo];
        CMz_Bridge = CMz_Bridge(KnInd);
        
        if opt.checks == 1       
                figure();
                semilogx(KnBridge,CD_Bridge,'-ok','LineWidth',2); hold on; grid on;
                stdPlot({'Kn','$C_D$'},{'$C_D$ Bridging'},{'Erf-based bridge','Sinusoidal Bridge'})
                xlim([1e-4 1e4])                             

                figure();
                semilogx(KnBridge,CL_Bridge,'-ok','LineWidth',2); hold on; grid on;
                stdPlot({'Kn','$C_L$'},{'$C_L$ Bridging'},{'Erf-based bridge','Sinusoidal Bridge'})
                xlim([1e-4 1e3])
                
                figure();
                semilogx(KnBridge,CS_Bridge,'-ok','LineWidth',2); hold on; grid on;
                stdPlot({'Kn','$C_S$'},{'$C_S$ Bridging'},{'Erf-based bridge','Sinusoidal Bridge'})
                xlim([1e-4 1e3])   

                figure();
                semilogx(KnBridge,CMx_Bridge,'-ok','LineWidth',2); hold on; grid on;
                stdPlot({'Kn','$C_{M,x}$'},{'$C_{M,x}$ Bridging'},{'Erf-based bridge','Sinusoidal Bridge'})
                xlim([1e-4 1e3])
                
                figure();
                semilogx(KnBridge,CMy_Bridge,'-ok','LineWidth',2); hold on; grid on;
                stdPlot({'Kn','$C_{M,y}$'},{'$C_{M,y}$ Bridging'},{'Erf-based bridge','Sinusoidal Bridge'})
                xlim([1e-4 1e3])
                
                figure();
                semilogx(KnBridge,CMz_Bridge,'-ok','LineWidth',2); hold on; grid on;
                stdPlot({'Kn','$C_{M,z}$'},{'$C_{M,z}$ Bridging'},{'Erf-based bridge','Sinusoidal Bridge'})
                xlim([1e-4 1e3])
                
                disp('transitional bridge check')
                keyboard
        end        
      
    end
    
%%     %%%%%%%%%%% TRANSLATIONAL HEATING %%%%%%%%%%%%%%%%%%%%
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.ThermoFlag == 1
%       preloading the bridges on built on for the referenced nose radius (Blunted cone 5deg, Mars Micro Probe, Mars Pathfinder, An averaged radius and Orion)
%         BridgeR = ppval(pchip_5Cone,Kn_trans_h);                              %% Blunted cone 5deg
        Rmodels(1) = []; % Removing the blunted cone
        BridgeR(1:length(Kn_trans_h),1) = ppval(pchip_Micro,Kn_trans_h);      %% Mars Micro Probe
%         Rmodels(2) = []; % Removing the mars microprobe
        BridgeR(1:length(Kn_trans_h),2) = ppval(pchip_MarsPath,Kn_trans_h);   %% Pathfinder
        BridgeR(1:length(Kn_trans_h),3) = ppval(pchip_MeanR,Kn_trans_h);      %% Average Rn
        BridgeR(1:length(Kn_trans_h),4) = ppval(pchip_Orion,Kn_trans_h);      %% Orion CEV
        BridgeR(BridgeR < 0) = 0;
        % building the global bridge for all the faces(and radius) and the requested Knudsen Number
        rN_bridge = rN;
        rN_bridge(rN_bridge > 5.3) = 5.3; % The maximum calibrated radius is 5.3m.
%         rN_bridge(rN_bridge < 0.0254) = 0.0254; % The minimum calibrated radius is 0.0254m. (Blunted Cone)
        rN_bridge(rN_bridge < 0.0875) = 0.0875; % The minimum calibrated radius is 0.0875m. (Mars Micro Probe)
        BridgeReq = pchip(Rmodels, BridgeR(:,:),rN_bridge)';
        
        
        % computing bridging Stc for each facet
        Stc_trans = repmat(Stc(:,end),1,length(Kn_trans_h)) + repmat(Stc_FMF(:,1)-Stc(:,end),1,length(Kn_trans_h)) .* BridgeReq;   
        Stc_Bridge = [Stc Stc_trans(:,2:end-1) Stc_FMF];
   
        if opt.checks == 1;
            [~,indToPlot] = sort(Stc(:,1));
            [~,indToPlotFMF] = sort(Stc_FMF(:,1));
            [~, RminInd] = min(rN);
            [~, RmaxInd] = max(rN);
            [~, thetaMaxInd] = max(Theta);

            figure();
            semilogx([Knc_h; Kn_trans_h; Knfmf],[Stc(opt.stagFaceIndvect(1),:) Stc_trans(opt.stagFaceIndvect(1),:) ...
                Stc_FMF(opt.stagFaceIndvect(1),:)],'-k','linewidth',3);  hold on; grid on;
            semilogx([Knc_h; Kn_trans_h; Knfmf],[Stc(thetaMaxInd,:) Stc_trans(thetaMaxInd,:) Stc_FMF(thetaMaxInd,:)],'--k','linewidth',3); 
            semilogx([Knc_h; Kn_trans_h; Knfmf],[Stc(RminInd,:) Stc_trans(RminInd,:) Stc_FMF(RminInd,:)],'-.k','linewidth',3);  
            semilogx([Knc_h; Kn_trans_h; Knfmf],[Stc(RmaxInd,:) Stc_trans(RmaxInd,:) Stc_FMF(RmaxInd,:)],':k','linewidth',3);
            stdPlot({'Kn','Stanton Number'},{'Stanton Number Bridging: radius-driven and normals-driven facets'},...
                {'Stagn. Point $\theta \simeq 0$','$\theta \simeq 90$','Min Radius','Max Radius'})
            xlim([1e-4 1e3]);
                        
            if FMFflag ==1 || bridgeFlag ==1
                figure(); hold on; grid on;
                trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
                trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Stc_FMF(:,1),'edgealpha',0.15); c = colorbar; set(gca, 'CLim', [0 1]) ;
                stdPlot({'x','y','z'},{'Fay-Riddel Local Radius-based Continuum Aero-thermal Heating'},{'off'},{['Kn = ',num2str(Knfmf(1))]})
            end
            if contFlag == 1 || bridgeFlag ==1;
                figure(); hold on; grid on;
                trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
                trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Stc(:,end),'edgealpha',0.15); c = colorbar; set(gca, 'CLim', [0 1]) ;
                stdPlot({'x','y','z'},{'Fay-Riddel Local Radius-based Continuum Aero-thermal Heating'},{'off'},{['Kn = ',num2str(Knc_h(end))]})
            end
            if bridgeFlag ==1
                figure(); hold on; grid on;
                trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.7 0.7 0.7]); hold on; axis('equal'); xlabel('X'); ylabel('Y'); zlabel('Z');
                trisurf(Fbfc,V(:,1),V(:,2),V(:,3),Stc_trans(:,2),'edgealpha',0.15); c = colorbar; set(gca, 'CLim', [0 1]) ;
                stdPlot({'x','y','z'},{'Fay-Riddel Local Radius-based Continuum Aero-thermal Heating'},{'off'},{['Kn = ',num2str(Kn_trans(2))]})
            end    
        figure(); hold on; grid on;
        plot(rN(indToPlot)/max(rN),'xb')        
        plot(Theta(indToPlot)/max(Theta),'ko')
        plot(Stc(indToPlot,1)/max(Stc(indToPlot,1)),'r','LineWidth',6)
        stdPlot({'Facet ID (sorted with increasing Stc)','Normalized Inputs'},{'Normalized plots: continuum regime facet Stc-Radius-$\theta$ correlation'}, ...
            {'Local Radius','Facet inclination $\theta$','F-R Stc'})
        
        figure(); hold on; grid on;
        plot(rN(indToPlotFMF)/max(rN),'bx')
        plot(Theta(indToPlotFMF)/max(Theta),'k','LineWidth',6)
        plot(Stc_FMF(indToPlotFMF,1)/max(Stc_FMF(indToPlotFMF,1)),'r','LineWidth',6)        
        stdPlot({'Facet ID (sorted with increasing Stc)','Normalized Inputs'},{'Normalized plots: FMF regime facet Stc-Radius-$\theta$ correlation'}, ...
            {'Local Radius','Facet Inclination $\theta$ ','FMF Heating'})
             keyboard
        end
    end %Thermo dynamics module ends
else
    if opt.AeroFlag == 1 && FMFflag == 1
        CDfmf = CDfm_memo';
        CLfmf =CLfm_memo';
        CSfmf = CSfm_memo';
        CMfmf = CMyfm_memo';     
    elseif opt.AeroFlag == 1 && contFlag == 1
        CDfmf = CDc_memo';
        CLfmf =CLc_memo';
        CSfmf = CSc_memo';
        CMfmf = CMyc_memo';     
    end
end

%% Simulation finished providing results depending on the used modules


disp(['simulation completed in: ',num2str(toc(tstart)),'s'])

IndOut = ismember(hInd,hIndReq);

if opt.AeroFlag == 1 && opt.ThermoFlag == 1
    if FMFflag == 1
        strath = {h(hIndReq), Kn(hIndReq) CDfm_memo(IndOut)' CLfm_memo(IndOut)' CSfm_memo(IndOut)' CMxfm_memo(IndOut)' ...
            CMyfm_memo(IndOut)' CMzfm_memo(IndOut)' StConst(:,IndOut) {F,V,N,preBFCind, N_verts} Stc_FMF(:,(IndOut))};
    elseif contFlag == 1
         strath = {h(hIndReq) Kn(hIndReq) CDc_memo(IndOut)' CLc_memo(IndOut)' CSc_memo(IndOut)' CMxc_memo(IndOut)' ...
             CMyc_memo(IndOut)' CMzc_memo(IndOut)' StConst(:,IndOut) {F,V,N,preBFCind, N_verts} Stc(:,(IndOut))};
    elseif FMFflag == 0 && contFlag == 0
        IndOut = ismember(hInd,hIndReq);
        strath = {h(hIndReq) Kn(hIndReq) CD_Bridge(IndOut)' CL_Bridge(IndOut)' CS_Bridge(IndOut)' CMx_Bridge(IndOut)' ...
            CMy_Bridge(IndOut)' CMz_Bridge(IndOut)' StConst(:,IndOut) {F,V,N,preBFCind, N_verts} Stc_Bridge(:,(IndOut))};
    end    
elseif opt.AeroFlag == 1 && opt.ThermoFlag == 0
%         disp('Aerothermodynamics has not been computed, the module was set to off')
    if FMFflag == 1
            strath = {h(hIndReq) Kn(hIndReq) CDfm_memo(IndOut)' CLfm_memo(IndOut)' CSfm_memo(IndOut)' CMxfm_memo(IndOut)' ...
                CMyfm_memo(IndOut)' CMzfm_memo(IndOut)'};
    elseif contFlag == 1
            strath = {h(hIndReq) Kn(hIndReq) CDc_memo(IndOut)' CLc_memo(IndOut)' CSc_memo(IndOut)' CMxc_memo(IndOut)' ...
                CMyc_memo(IndOut)' CMzc_memo(IndOut)'};
    elseif FMFflag == 0 && contFlag == 0
            strath = {h(hIndReq) Kn(hIndReq) CD_Bridge(IndOut)' CL_Bridge(IndOut)' CS_Bridge(IndOut)' CMx_Bridge(IndOut)' ...
                CMy_Bridge(IndOut)' CMz_Bridge(IndOut)'};
    end    
elseif opt.AeroFlag == 0 && opt.ThermoFlag == 1
%     disp('Aerodynamics have not been computed, the module was set to off')
    strath = {h(hIndReq) Kn(hIndReq) 0 CH_stagReq Kn(hIndReq)};
end


end
