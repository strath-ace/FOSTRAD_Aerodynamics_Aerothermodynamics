%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOSTRAD: The Free Open Source Tool for Re-entry of Asteroids and Debris %%
%%%%%%%%  Copyright (C) 2021 University of Strathclyde and Authors %%%%%%%%%%
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

function [varargout] = refineSTL_tri_v2_2(F, V,checkType, MaxPar,MAXFACES,NmaxIter,shadowingRefFlag,figFlag)
% refineSTL_tri_v2_1(F, V,checkType, MaxPar,MAXFACES,NmaxIter,shadowingRefFlag,figFlag)
% Refines a given STL geometry according to maximum triangles area or side lengths
%
% outputs: argout = {F V N};
% Inputs: 
% - F           set of connections to define the F [Nx3] where N is the number of F
% - V           set of points coordinates in the form of [Mx3] [x y z] where M is the number of vertexes
% - checktype   string to define whether to use the <'area'> or the <'length'> or <'both'> to refine 
%               by both: the area and the length
% - MaxPar      Parameter to define the maximum area or length to be refined. In case the STL geometry
%               is being refined with <'both'> MaxPar is a a vector defined as [MaxArea MaxLength]
% - MAXFACES    By default the mesh refinement is stopped as it reaches 20'000 F. Use this parameter
%               if you want a higher refinement. WARNING: too many F could saturate the RAM
%               The parameter is very useful also when the refinement must not exceed an arbitrary F
%               Number, e.g.: when you want to limit the computational expense
% - NmaxIter    Number of iterations performed to refine the mesh. At every iteration it refines a 
%               (NmaxFaces-NstartingFaces)/iter faces, giving priority to the faces with higher area
%               or length. In case <'both'> refinement is active, the algorithm will prioritize the area.
% - shadowingRefFlag       When coupled with an occlusion culling algorithm, the refinement can be applied 
%                          to the shadowed-non-shadowed bordering facets. In order to use this feature
%                          the user must specify the [shadowingRefFlag] = {'on',[refIndex]} where
%                          [refIndex] is the matrix containing the indexes of the faces to be refined
%                          based on the [F] indexing reference. The code will adaptively refine the facets
%                          given as input with the normal procedure.
% - figFlag     Just enter any value to activate the fig, remember to turn off the shadowing refinement.
%               e.g.: refineSTL_tri_v2_1(F, V,checkType, MaxPar,MAXFACES,NmaxIter,{'off'},['ON'])
if exist('shadowingRefFlag','var')
    if strcmpi(shadowingRefFlag{1},'on')
        if length(shadowingRefFlag) < 2
            error('You must provide the indexes to be refined on the F reference. shadowingRefFlag = {''on'', indexToRef}')
        else
            indexToRefReq = shadowingRefFlag{2};
        end
        shadowingRefFlag = 1;
    else
        shadowingRefFlag = 0;
    end
else
    shadowingRefFlag = 0;
end

combRef.boh(1,1:2) = [1 2]; % 1 3 2 3];
combRef.boh(2,1:2) = [1 3];
combRef.boh(3,1:2) = [2 3];
if exist('MAXFACES','var')
    maxFaces = MAXFACES;
else
    maxFaces = 2e4;
end

NstartFaces = length(F);
NfacesPerIter = ceil((maxFaces-NstartFaces)/NmaxIter/2);
counter = 0;
tstart = tic;
% Amult = 5;
    %% Area Refinement
    
if strcmpi(checkType,'area')
    AreaMax = MaxPar;
    nVertOld = size(V,1); 
    Nfaces = length(F);
    
    % find the indexes to refine by Area
    area_Tri = Area(F, V);
    IndToRef = find(area_Tri > AreaMax);
    [~, indSortA] = sort(area_Tri(IndToRef),'descend');
    IndToRef = IndToRef(indSortA);

    if shadowingRefFlag ==1
        [IndToRef,~,IndReqToUpdate] = intersect(IndToRef, indexToRefReq,'stable');
        if length(IndToRef) >=  NfacesPerIter
            IndToRef = IndToRef(1:NfacesPerIter);
        end
        IndReqToUpdate = IndReqToUpdate(1:length(IndToRef));
        IndReqNotToUpdate = setdiff(1:length(indexToRefReq),IndReqToUpdate)';
        IndReqNotToUpdate = indexToRefReq(IndReqNotToUpdate);
    else
        if length(IndToRef) >=  NfacesPerIter
            IndToRef = IndToRef(1:NfacesPerIter);
        end
        IndReqNotToUpdate = [];
    end    
        
    while counter <= NmaxIter
        clear LogicMat facesNew
        
        e0=V(F(IndToRef,3),:)-V(F(IndToRef,2),:);
        e1=V(F(IndToRef,1),:)-V(F(IndToRef,3),:);
        N=cross(e1,e0);
        N=normr(N); 
        N=repmat(N,2,1);
        
        
        IndnotRef = setdiff(1:length(F),IndToRef)';
        Nind = size(IndToRef,1);
        CombTemp = repmat(combRef,Nind,1,1);
        Dxyz12 = V(F(IndToRef,1),:)-V(F(IndToRef,2),:);
        Dxyz13 = V(F(IndToRef,1),:)-V(F(IndToRef,3),:);
        Dxyz23 = V(F(IndToRef,2),:)-V(F(IndToRef,3),:);
        L12 = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
        L13 = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
        L23 = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
        IndCombID = sum((repmat(max([L12 L13 L23],[],2),1,3) == [L12 L13 L23]).*[ones(Nind,1) ones(Nind,1)*2 ones(Nind,1)*3],2);
        IndCombID(IndCombID==5) = 3;
        IndCombID(IndCombID==4) = 2;
        Ind2sel = [combRef.boh(IndCombID,1) combRef.boh(IndCombID,2)];
        LogicMat(Ind2sel(:,1)==1,1) = 1;
        LogicMat(Ind2sel(:,1)==2 | Ind2sel(:,2)==2,2) = 1;
        LogicMat(Ind2sel(:,2)==3,3) = 1;
        LogicMat = logical(LogicMat)';
        F_new = reshape(F(IndToRef,:),[],3)';
        newVertInd = reshape(F_new(LogicMat),2,[]);
        newVert = [V(newVertInd(1,:),:)  V(newVertInd(2,:),:)];
        newVert = round((newVert(:,1:3)+newVert(:,4:6))/2,9);
        MissInd = 6-sum(Ind2sel,2);
%         for i = 1:Nind
%            facesNew(i*2-1:i*2,:) = [F(IndToRef(i),MissInd(i)) F(IndToRef(i),Ind2sel(i,1)) i+nVertOld; F(IndToRef(i),MissInd(i)) F(IndToRef(i),Ind2sel(i,2)) i+nVertOld];
%         end
        facesNew = [F(IndToRef+(MissInd(:)-1)*Nfaces) F(IndToRef+(Ind2sel(:,1)-1)*Nfaces) [1:Nind]'+nVertOld
                    F(IndToRef+(MissInd(:)-1)*Nfaces) F(IndToRef+(Ind2sel(:,2)-1)*Nfaces) [1:Nind]'+nVertOld];
        NfacesNew = size(facesNew,1);
        FtempNotToUpdate = F(IndReqNotToUpdate,:);
        V = [V; newVert]; 
        
        % correcting the winding direction
        e0=V(facesNew(:,3),:)-V(facesNew(:,2),:);
        e1=V(facesNew(:,1),:)-V(facesNew(:,3),:);
        Nnew=cross(e1,e0);
        Nnew=normr(Nnew);   
        CorrectN = dot(Nnew,N,2);
        facesNewTemp = facesNew;
        facesNew(CorrectN<0,2) = facesNewTemp(CorrectN<0,3);
        facesNew(CorrectN<0,3) = facesNewTemp(CorrectN<0,2);        
        
        
        F = [F(IndnotRef,:); facesNew];         
        nVertOld = size(V,1);  
        Nfaces = length(F);
        
        
        if shadowingRefFlag == 1         
            [~,~,indexToRefReq] = intersect(FtempNotToUpdate,F,'rows','stable');
            indexToRefReq(length(indexToRefReq)+1:length(indexToRefReq)+NfacesNew) = Nfaces-NfacesNew+1:Nfaces;
        end
        
        
        % find the indexes to refine by Area
        area_Tri = Area(F, V);
        IndToRef = find(area_Tri > AreaMax);
        [~, indSortA] = sort(area_Tri(IndToRef),'descend');
        IndToRef = IndToRef(indSortA);
        
        if shadowingRefFlag ==1
            [IndToRef,~,IndReqToUpdate] = intersect(IndToRef, indexToRefReq,'stable');
            if length(IndToRef) >=  NfacesPerIter
                IndToRef = IndToRef(1:NfacesPerIter);
            end
            IndReqToUpdate = IndReqToUpdate(1:length(IndToRef));
            IndReqNotToUpdate = setdiff(1:length(indexToRefReq),IndReqToUpdate)';
            IndReqNotToUpdate = indexToRefReq(IndReqNotToUpdate);
        else
            if length(IndToRef) >=  NfacesPerIter
                IndToRef = IndToRef(1:NfacesPerIter);
            end
            IndReqNotToUpdate = [];
        end        

        counter = counter +1;
        NfacesPerIter = ceil((maxFaces-Nfaces)/(NmaxIter-counter)/2);
        if isempty(IndToRef)
            disp(['Maximum area refinement reached - no more faces to refined - decrease the max area for a better refinement'])
            break
        end            
        if length(F) > Nfaces
            disp(['Number of Maximum Facets [',num2str(maxFaces),'] reached. Increase the Maximum number of F to increase the refinement'])
            break
        end

    end
    disp(['Triangular Mesh Refinement finished in ',num2str(counter-1),' iterations and ',num2str(toc(tstart)),'s'])
    if exist('figFlag','var')    
        figure();
        trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; xlabel('X'); ylabel('Y'); zlabel('Z'); axis('equal');
    end

    
    %% Length Refinement
    
elseif strcmpi(checkType, 'length')
    nVertOld = size(V,1);
    Nfaces = length(F);
    Lmax = MaxPar(1);
    
    % find the indexes to refine by length
    Dxyz12 = V(F(:,1),:)-V(F(:,2),:);
    Dxyz13 = V(F(:,1),:)-V(F(:,3),:);
    Dxyz23 = V(F(:,2),:)-V(F(:,3),:);
    L(1:length(Dxyz12),1) = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
    L(:,2) = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
    L(:,3) = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
    LmaxV = max(L,[],2);
    IndToRef = find(max(L,[],2) > Lmax);
    [~, indSort] = sort(LmaxV(IndToRef),'descend');
    IndToRef = IndToRef(indSort);
    
    if shadowingRefFlag ==1
        [IndToRef,~,IndReqToUpdate] = intersect(IndToRef, indexToRefReq,'stable');
        if length(IndToRef) >=  NfacesPerIter
            IndToRef = IndToRef(1:NfacesPerIter);
            %             trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
            %             trisurf(F(IndToRef,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]); hold on; grid on; axis('equal');
        end
        IndReqToUpdate = IndReqToUpdate(1:length(IndToRef));
        IndReqNotToUpdate = setdiff(1:length(indexToRefReq),IndReqToUpdate)';
        IndReqNotToUpdate = indexToRefReq(IndReqNotToUpdate);
    else
        if length(IndToRef) >=  NfacesPerIter
            IndToRef = IndToRef(1:NfacesPerIter);
        end
        IndReqNotToUpdate = [];
    end
 
    
    while counter <= NmaxIter
        clear LogicMat facesNew

        e0=V(F(IndToRef,3),:)-V(F(IndToRef,2),:);
        e1=V(F(IndToRef,1),:)-V(F(IndToRef,3),:);
        N=cross(e1,e0);
        N=normr(N); 
        N=repmat(N,2,1);
        
        IndnotRef = setdiff(1:length(F),IndToRef)';
        Nind = size(IndToRef,1);
        CombTemp = repmat(combRef,Nind,1,1);
        Dxyz12 = V(F(IndToRef,1),:)-V(F(IndToRef,2),:);
        Dxyz13 = V(F(IndToRef,1),:)-V(F(IndToRef,3),:);
        Dxyz23 = V(F(IndToRef,2),:)-V(F(IndToRef,3),:);
        L12 = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
        L13 = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
        L23 = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
        IndCombID = sum((repmat(max([L12 L13 L23],[],2),1,3) == [L12 L13 L23]).*[ones(Nind,1) ones(Nind,1)*2 ones(Nind,1)*3],2);
        IndCombID(IndCombID==5) = 3;
        IndCombID(IndCombID==4) = 2;
        Ind2sel = [combRef.boh(IndCombID,1) combRef.boh(IndCombID,2)];
        LogicMat(Ind2sel(:,1)==1,1) = 1;
        LogicMat(Ind2sel(:,1)==2 | Ind2sel(:,2)==2,2) = 1;
        LogicMat(Ind2sel(:,2)==3,3) = 1;
        LogicMat = logical(LogicMat)';
        F_new = reshape(F(IndToRef,:),[],3)';
        newVertInd = reshape(F_new(LogicMat),2,[]);
        newVert = [V(newVertInd(1,:),:)  V(newVertInd(2,:),:)];
        newVert = round((newVert(:,1:3)+newVert(:,4:6))/2,9);
        MissInd = 6-sum(Ind2sel,2);
%         for i = 1:Nind
%            facesNew(i*2-1:i*2,:) = [F(IndToRef(i),MissInd(i)) F(IndToRef(i),Ind2sel(i,1)) i+nVertOld; F(IndToRef(i),MissInd(i)) F(IndToRef(i),Ind2sel(i,2)) i+nVertOld];
%         end
        facesNew = [F(IndToRef+(MissInd(:)-1)*Nfaces) F(IndToRef+(Ind2sel(:,1)-1)*Nfaces) [1:Nind]'+nVertOld
                    F(IndToRef+(MissInd(:)-1)*Nfaces) F(IndToRef+(Ind2sel(:,2)-1)*Nfaces) [1:Nind]'+nVertOld];
        NfacesNew = size(facesNew,1);
        FtempNotToUpdate = F(IndReqNotToUpdate,:);
        V = [V; newVert]; 
        
        % correcting the winding direction
        e0=V(facesNew(:,3),:)-V(facesNew(:,2),:);
        e1=V(facesNew(:,1),:)-V(facesNew(:,3),:);
        Nnew=cross(e1,e0);
        Nnew=normr(Nnew);   
        CorrectN = dot(Nnew,N,2);
        facesNewTemp = facesNew;
        facesNew(CorrectN<0,2) = facesNewTemp(CorrectN<0,3);
        facesNew(CorrectN<0,3) = facesNewTemp(CorrectN<0,2);
        
        
        F = [F(IndnotRef,:); facesNew];       
        nVertOld = size(V,1);
        Nfaces = length(F);
        
        if shadowingRefFlag == 1         
            [~,~,indexToRefReq] = intersect(FtempNotToUpdate,F,'rows','stable');
            indexToRefReq(length(indexToRefReq)+1:length(indexToRefReq)+NfacesNew) = Nfaces-NfacesNew+1:Nfaces;
        end
        
        % find the indexes to refine by length
        Dxyz12 = V(F(:,1),:)-V(F(:,2),:);
        Dxyz13 = V(F(:,1),:)-V(F(:,3),:);
        Dxyz23 = V(F(:,2),:)-V(F(:,3),:);
        L(1:length(Dxyz12),1) = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
        L(:,2) = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
        L(:,3) = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
        LmaxV = max(L,[],2);
        IndToRef = find(max(L,[],2) > Lmax);
        [~, indSort] = sort(LmaxV(IndToRef),'descend');
        IndToRef = IndToRef(indSort);
        if shadowingRefFlag ==1
            [IndToRef,~,IndReqToUpdate] = intersect(IndToRef, indexToRefReq,'stable');
            if length(IndToRef) >=  NfacesPerIter
                IndToRef = IndToRef(1:NfacesPerIter);
                %             trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
                %             trisurf(F(IndToRef,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]); hold on; grid on; axis('equal');
            end
            IndReqToUpdate = IndReqToUpdate(1:length(IndToRef));
            IndReqNotToUpdate = setdiff(1:length(indexToRefReq),IndReqToUpdate)';
            IndReqNotToUpdate = indexToRefReq(IndReqNotToUpdate);
        else
            if length(IndToRef) >=  NfacesPerIter
                IndToRef = IndToRef(1:NfacesPerIter);
            end
            IndReqNotToUpdate = [];
        end
        

        counter = counter +1;
        NfacesPerIter = ceil((maxFaces-Nfaces)/(NmaxIter-counter)/2);
        if isempty(IndToRef)
            disp(['Maximum side length refinement reached - no more faces to refined - decrease the max length for a better refinement'])
            break
        end       
        if Nfaces > maxFaces
            disp(['Number of Maximum Facets [',num2str(maxFaces),'] reached. Increase the Maximum number of F to increase the refinement'])
            break
        end        
        
    end
    disp(['Triangular Mesh Refinement finished in ',num2str(counter-1),' iterations and ',num2str(toc(tstart)),'s'])
    if exist('figFlag','var')
        figure();
        trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on;  xlabel('X'); ylabel('Y'); zlabel('Z');axis('equal');
    end
    
    %% Area and Length Refinement
elseif strcmpi(checkType, 'both') && length(MaxPar) == 2
    nVertOld = size(V,1);
    Nfaces = length(F);
    AreaMax = MaxPar(1);
    Lmax = MaxPar(2);
    
    % find the indexes to refine by area
    area_Tri = Area(F, V);
    IndToRefA = find(area_Tri > AreaMax);
    [~, indSortA] = sort(area_Tri(IndToRefA),'descend');
    IndToRefA = IndToRefA(indSortA);
    
    

    % find the indexes to refine by length    
    Dxyz12 = V(F(:,1),:)-V(F(:,2),:);
    Dxyz13 = V(F(:,1),:)-V(F(:,3),:);
    Dxyz23 = V(F(:,2),:)-V(F(:,3),:);
    L(1:length(Dxyz12),1) = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
    L(:,2) = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
    L(:,3) = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
    LmaxV = max(L,[],2);
    IndToRefL = find(LmaxV > Lmax);
    [~, indSortL] = sort(LmaxV(IndToRefL),'descend');
    IndToRefL = IndToRefL(indSortL);
    

    
    IndToRef = unique([IndToRefA; IndToRefL],'stable');
    if shadowingRefFlag ==1
        [IndToRef,~,IndReqToUpdate] = intersect(IndToRef, indexToRefReq,'stable');
        if length(IndToRef) >=  NfacesPerIter
            IndToRef = IndToRef(1:NfacesPerIter);
            %             trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
            %             trisurf(F(IndToRef,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]); hold on; grid on; axis('equal');
        end
        IndReqToUpdate = IndReqToUpdate(1:length(IndToRef));
        IndReqNotToUpdate = setdiff(1:length(indexToRefReq),IndReqToUpdate)';
        IndReqNotToUpdate = indexToRefReq(IndReqNotToUpdate);
    else
        if length(IndToRef) >=  NfacesPerIter
            IndToRef = IndToRef(1:NfacesPerIter);
        end
        IndReqNotToUpdate = [];
    end
    
    while counter <= NmaxIter
        clear LogicMat facesNew
        
        e0=V(F(IndToRef,3),:)-V(F(IndToRef,2),:);
        e1=V(F(IndToRef,1),:)-V(F(IndToRef,3),:);
        N=cross(e1,e0);
        N=normr(N); 
        N=repmat(N,2,1);
        
        IndnotRef = setdiff(1:length(F),IndToRef)';
        Nind = size(IndToRef,1);
        CombTemp = repmat(combRef,Nind,1,1);
        Dxyz12 = V(F(IndToRef,1),:)-V(F(IndToRef,2),:);
        Dxyz13 = V(F(IndToRef,1),:)-V(F(IndToRef,3),:);
        Dxyz23 = V(F(IndToRef,2),:)-V(F(IndToRef,3),:);
        L12 = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
        L13 = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
        L23 = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
        IndCombID = sum((repmat(max([L12 L13 L23],[],2),1,3) == [L12 L13 L23]).*[ones(Nind,1) ones(Nind,1)*2 ones(Nind,1)*3],2);
        IndCombID(IndCombID==5) = 3;
        IndCombID(IndCombID==4) = 2;
        Ind2sel = [combRef.boh(IndCombID,1) combRef.boh(IndCombID,2)];
        LogicMat(Ind2sel(:,1)==1,1) = 1;
        LogicMat(Ind2sel(:,1)==2 | Ind2sel(:,2)==2,2) = 1;
        LogicMat(Ind2sel(:,2)==3,3) = 1;
        LogicMat = logical(LogicMat)';
        F_new = reshape(F(IndToRef,:),[],3)';
        newVertInd = reshape(F_new(LogicMat),2,[]);
        newVert = [V(newVertInd(1,:),:)  V(newVertInd(2,:),:)];
        newVert = round((newVert(:,1:3)+newVert(:,4:6))/2,9);
        MissInd = 6-sum(Ind2sel,2);
%         for i = 1:Nind
%            facesNew(i*2-1:i*2,:) = [F(IndToRef(i),MissInd(i)) F(IndToRef(i),Ind2sel(i,1)) i+nVertOld; F(IndToRef(i),MissInd(i)) F(IndToRef(i),Ind2sel(i,2)) i+nVertOld];
%         end
        facesNew = [F(IndToRef+(MissInd(:)-1)*Nfaces) F(IndToRef+(Ind2sel(:,1)-1)*Nfaces) [1:Nind]'+nVertOld
                    F(IndToRef+(MissInd(:)-1)*Nfaces) F(IndToRef+(Ind2sel(:,2)-1)*Nfaces) [1:Nind]'+nVertOld];    
        NfacesNew = size(facesNew,1);
        FtempNotToUpdate = F(IndReqNotToUpdate,:);
        V = [V; newVert];
        
        % correcting the winding direction
        e0=V(facesNew(:,3),:)-V(facesNew(:,2),:);
        e1=V(facesNew(:,1),:)-V(facesNew(:,3),:);
        Nnew=cross(e1,e0);
        Nnew=normr(Nnew);   
        CorrectN = dot(Nnew,N,2);
        facesNewTemp = facesNew;
        facesNew(CorrectN<0,2) = facesNewTemp(CorrectN<0,3);
        facesNew(CorrectN<0,3) = facesNewTemp(CorrectN<0,2);
        
        F = [F(IndnotRef,:); facesNew];        
        Nfaces = size(F,1);
        nVertOld = size(V,1);
        
        
        if shadowingRefFlag == 1         
            [~,~,indexToRefReq] = intersect(FtempNotToUpdate,F,'rows','stable');
            indexToRefReq(length(indexToRefReq)+1:length(indexToRefReq)+NfacesNew) = Nfaces-NfacesNew+1:Nfaces;
%             trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
%             trisurf(F(indexToRefReq,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0 0 1]); hold on; grid on; axis('equal');
        end

    % find the indexes to refine by area
    area_Tri = Area(F, V);
    IndToRefA = find(area_Tri > AreaMax);
    [~, indSortA] = sort(area_Tri(IndToRefA),'descend');

    % find the indexes to refine by length    
    Dxyz12 = V(F(:,1),:)-V(F(:,2),:);
    Dxyz13 = V(F(:,1),:)-V(F(:,3),:);
    Dxyz23 = V(F(:,2),:)-V(F(:,3),:);
    L(1:length(Dxyz12),1) = sqrt(Dxyz12(:,1).^2+Dxyz12(:,2).^2+Dxyz12(:,3).^2);
    L(:,2) = sqrt(Dxyz13(:,1).^2+Dxyz13(:,2).^2+Dxyz13(:,3).^2);
    L(:,3) = sqrt(Dxyz23(:,1).^2+Dxyz23(:,2).^2+Dxyz23(:,3).^2);
    LmaxV = max(L,[],2);
    IndToRefL = find(LmaxV > Lmax);
    [~, indSortL] = sort(LmaxV(IndToRefL),'descend');
    IndToRefL = IndToRefL(indSortL);
    
    
    IndToRef = unique([IndToRefA; IndToRefL],'stable');
        if shadowingRefFlag ==1
            [IndToRef,~,IndReqToUpdate] = intersect(IndToRef, indexToRefReq,'stable');
            if length(IndToRef) >=  NfacesPerIter
                IndToRef = IndToRef(1:NfacesPerIter);
                %             trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on; grid on; axis('equal');
                %             trisurf(F(IndToRef,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[1 0 0]); hold on; grid on; axis('equal');
            end
            IndReqToUpdate = IndReqToUpdate(1:length(IndToRef));
            IndReqNotToUpdate = setdiff(1:length(indexToRefReq),IndReqToUpdate)';
            IndReqNotToUpdate = indexToRefReq(IndReqNotToUpdate);
        else
            if length(IndToRef) >=  NfacesPerIter
                IndToRef = IndToRef(1:NfacesPerIter);
            end
            IndReqNotToUpdate = [];
        end
        
%         tr = triangulation(F,V);
%         faceBaricenter = incenter(tr);
%         e0=V(F(:,2),:)-V(F(:,3),:);
%         e1=V(F(:,1),:)-V(F(:,2),:);
%         N=cross(e1,e0);        
%         if shadowingRefFlag == 1
%             Az = linspace(-135,-45,NmaxIter+1);
%             f1 = figure('Position', [100, 100, 800, 600],'visible','off');
% %             f1 = figure('Position', [100, 100, 800, 600]);
%             trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0.65 0.65 0.65]); hold on;  xlabel('X'); ylabel('Y'); zlabel('Z'); axis('equal');
%             trisurf(F(indexToRefReq,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.25,'FaceColor',[0 0 1]); hold on;  xlabel('X'); ylabel('Y'); zlabel('Z'); axis('equal');        
% 
%     %         fq1 = quiver3(faceBaricenter(:,1),faceBaricenter(:,2),faceBaricenter(:,3), N(:,1), N(:,2), N(:,3),3,'r');
%     %         view(Az(counter+1),25);
%             view(-110,25)
%             if exist('AdaptiveFrames.mat')
%                 if counter == 0;
%                     load('AdaptiveFrames.mat')
%                 end
%                 Frames(length(Frames)+1)= getframe(f1);
%             else
%                 Frames(1+counter)= getframe(f1);;
%             end
%             close(f1)
%         end
        
        counter = counter +1;
        NfacesPerIter = ceil((maxFaces-Nfaces)/(NmaxIter-counter)/2);
        if isempty(IndToRef)
            disp(['Maximum side length or area refinement reached - no more faces to refined - decrease the max area/length for a better refinement'])
            break
        end             
        if Nfaces > maxFaces
            disp(['Number of Maximum Facets [',num2str(maxFaces),'] reached. Increase the Maximum number of F to increase the refinement'])
            break
        end        
        
    end
    disp(['Triangular Mesh Refinement finished in ',num2str(counter-1),' iterations and ',num2str(toc(tstart)),'s'])
%     if shadowingRefFlag == 1
%         save('AdaptiveFrames.mat','Frames')
%     end
    if exist('figFlag','var')
        figure();
        trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on;  xlabel('X'); ylabel('Y'); zlabel('Z'); axis('equal');
    end

elseif strcmpi(checkType, 'both') && length(MaxPar) == 1
    error('You must enter two different value in the form of MaxPar = [AreaMax Lmax]')
else
    error('The Check type you entered is invalid, select between checkType = <''area''> or <''length''> or <''both''>' )
end

e0=V(F(:,3),:)-V(F(:,2),:);
e1=V(F(:,1),:)-V(F(:,3),:);
N=cross(e1,e0);
N=normr(N);

if nargout == 3;
    varargout = {F V N};
else
    varargout = {F V};
end

end

