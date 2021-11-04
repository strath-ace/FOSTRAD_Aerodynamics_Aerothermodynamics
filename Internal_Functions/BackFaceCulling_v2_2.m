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

function varargout = BackFaceCulling_v2_2(AofA, SS, Roll, F, V, N, Npix, whiteFlag, BFCflag, shellFlag)
% profile on
% Computing the backface culling
% [bfcFVN{1}, bfcFVN{2}, bfcFVN{3}] = BackFaceCulling_v1(AofA, SS, F, V, N, Npix, whiteFlag, BFCflag);
% - AofA = angle of attack (pitch)
% - SS = side slip angle (yaw)
% - F = matrix containing the links for the faces
% - V = vertixes of the connected points [X, Y, Z]
% - N = normals
% - if you want to compute only the wetted facets according to the wind direction
% - number of pixel for the images resolution, increases the accuracy and computational time
% 
% % plotting the backface culling
% figure('Position', [100, 100, 1049, 895]);
% trisurf(bfcFVN{1}, bfcFVN{2}(:,1),bfcFVN{2}(:,2),bfcFVN{2}(:,3),'EdgeAlpha',0.15); view([90-SS, AofA])
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('BackFaceCulling Resulting mesh [m]')

                                            %% BACKFACE CULLING ALGORITHM %%
tstart = tic;
% close all

% cleaning and re-ordening the mesh
% [V, F] = reOrderMesh(V, F, 3);                                        
% 
Edge = 0.05;

% % % filtering survi
% % F = F(surv_IDs,:);
% % N = N(surv_IDs,:);

% detect wetted triangles according to their normals;                                        
ROT = Rotation_v1(SS,AofA,Roll);
if BFCflag == 1 && shellFlag == 0
    WindDir = [ROT*[1 0 0]']';
    WettedInd = dot(repmat(WindDir,length(N),1),N,2);
    WettedInd = find(WettedInd < Edge);
    N = N(WettedInd,:);
    F = F(WettedInd,:);
%     trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]);
elseif shellFlag == 1
    WindDir = [ROT*[1 0 0]']';
    WettedInd = dot(repmat(WindDir,length(N),1),N,2);
    WettedInd = find(WettedInd < Edge);
    InnerFaces = setdiff((1:length(F))',WettedInd);   
%     trisurf(F(WettedInd,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0.65 0.65]); hold on;
%     trisurf(F(InnerFaces,:),V(:,1),V(:,2),V(:,3),'edgealpha',0.15,'FaceColor',[0.65 0 0]);
end
% 
%     faceBaricenter = sum(reshape(V(F(:,1:3),1),[],3),2)/3;
%     faceBaricenter(:,2) = sum(reshape(V(F(:,1:3),2),[],3),2)/3;
%     faceBaricenter(:,3) = sum(reshape(V(F(:,1:3),3),[],3),2)/3;
%     fq1 = quiver3(faceBaricenter(:,1),faceBaricenter(:,2),faceBaricenter(:,3), N(:,1), N(:,2), N(:,3),3,'r'); 


% % directionX = [1 0 0];
% % directionY = [0 1 0];
% % directionZ = [0 0 1];
% % % G = randi([0, 255],length(F),3);

Nfaces = length(F);
steps = ceil(Nfaces^(1/3));
% 
% Div = [25 20 15 10 5 3 1];
% BaseDiv = floor(32./Div);
% Ncomb = BaseDiv.^3;
% indToSel = find(Ncomb > Nfaces);
if whiteFlag == 1 
    comb = floor(linspace(1,steps,steps+1));
    comb(comb==240) = 241;
    [G1, G2, G3] = ndgrid(comb, comb, comb);
    G = [reshape(G1,[],1) reshape(G2,[],1) reshape(G3,[],1)];
    G = G(1:Nfaces,:);
else
    G = randi([0, 255],length(F),3);
end

cycle = 0;

[Gu, ind] = unique(G, 'rows');
duplicate_ind = setdiff([1:size(G, 1)]', ind);
Gblack = find(sum(G,2)==0);
if isempty(Gblack) == 0
    Gsub = randi([0, 255],length(Gblack),3);
    for i = 1:length(Gblack)
        G(Gblack(i)) = Gsub(i);
    end
end


while length(Gu) ~= length(G) && ~isempty(duplicate_ind) && cycle < 100 %% fix the issue here!!!
    G(duplicate_ind,:) = randi([0, 255],length(duplicate_ind),3);
    [Gu, ind] = unique(G, 'rows');
    duplicate_ind = setdiff([1:size(G, 1)]', ind);
    cycle = cycle +1;
end

G = round(G/255,4);
ColorID = uint8(G*255);

f1 = figure('Position', [0, 0, Npix, Npix],'visible','off');
% f1 = figure('Position', [100, 100, Npix, Npix]);
% axis('equal')
Vt = (ROT'*V')';
if whiteFlag == 1
    p=patch('Faces',F,'Vertices',Vt,'FaceVertexCdata',G,'FaceColor','flat','CDataMapping','direct','LineWidth',0.01,'EdgeColor',[1 1 1]);
else
    p=patch('Faces',F,'Vertices',Vt,'FaceVertexCdata',G,'FaceColor','flat','LineStyle','none','CDataMapping','direct');
end
xlabel('X')
ylabel('Y')
zlabel('Z')
% axis('equal')
hold on
view([-90,0,0])
axis off



imageData = getframe(f1);
close all
imageData = imageData.cdata;
% imshow(imageData)
imageData = reshape(imageData,[],1,3);
imageData = [imageData(:,1,1), imageData(:,1,2), imageData(:,1,3)];
imageData(imageData(:,1) == 240 & imageData(:,2) == 240 & imageData(:,3) == 240,:) = [];
imageData = unique(imageData, 'rows');
[~, ~, IB] = intersect(imageData, ColorID,'rows');



% Assigning the updated faces and normals
if shellFlag == 0
    F = F(IB,:);
    N = N(IB,:);
else
    [~, ~, InvertedFacesIDs] = intersect(InnerFaces,IB);
    SignNorm = ones(length(IB),1);
    SignNorm(InvertedFacesIDs) = -1;
    F = F(IB,:);
    N = N(IB,:).*SignNorm;
end
% figure();
% trisurf(F,V(:,1),V(:,2),V(:,3),'EdgeAlpha',0.15)
% axis('equal')
tcomp = toc(tstart);
% disp(['BackFaceCulling Shape initialised in: ', num2str(tcomp),'s'])

if nargout == 3;
    varargout = {F, V, N};
elseif nargout == 4;
    varargout = {F, V, N, tcomp};
end

% profile viewer

end

