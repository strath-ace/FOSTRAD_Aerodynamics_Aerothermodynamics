function [VertexNormals,Avertex,Acorner,up,vp, avEdge]=CalcVertexNormals(FV,N)
%% Summary
%Author: Itzik Ben Shabat
%Last Update: July 2014

%summary: CalcVertexNormals calculates the normals and voronoi areas at each vertex
%INPUT:
%FV - triangle mesh in face vertex structure
%N - face normals
%OUTPUT -
%VertexNormals - [Nv X 3] matrix of normals at each vertex
%Avertex - [NvX1] voronoi area at each vertex
%Acorner - [NfX3] slice of the voronoi area at each face corner
%% License
% Copyright (c) 2016, Itzik Ben Shabat
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of Technion nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%% Code
% disp('Calculating vertex normals... Please wait');


% Get all edge vectors
e0=FV.vertices(FV.faces(:,3),:)-FV.vertices(FV.faces(:,2),:);
e1=FV.vertices(FV.faces(:,1),:)-FV.vertices(FV.faces(:,3),:);
e2=FV.vertices(FV.faces(:,2),:)-FV.vertices(FV.faces(:,1),:);
% Normalize edge vectors
e0_norm=normr(e0);
e1_norm=normr(e1);
e2_norm=normr(e2);

%normalization procedure
%calculate face Area
%edge lengths
de0=sqrt(e0(:,1).^2+e0(:,2).^2+e0(:,3).^2);
de1=sqrt(e1(:,1).^2+e1(:,2).^2+e1(:,3).^2);
de2=sqrt(e2(:,1).^2+e2(:,2).^2+e2(:,3).^2);

% de0(de0<=eps(1e9)) = eps(1e9);
% de1(de1<=eps(1e9)) = eps(1e9);
% de2(de2<=eps(1e9)) = eps(1e9);
% 
De0 = mean(de0(de0>1e-6));
De1 = mean(de1(de1>1e-6));
De2 = mean(de2(de2>1e-6));


avEdge = mean([De0,De1,De2]);


l2=[de0.^2 de1.^2 de2.^2];

%using ew to calulate the cot of the angles for the voronoi area
%calculation. ew is the triangle barycenter, I later check if its inside or
%outide the triangle
ew=[l2(:,1).*(l2(:,2)+l2(:,3)-l2(:,1)) l2(:,2).*(l2(:,3)+l2(:,1)-l2(:,2)) l2(:,3).*(l2(:,1)+l2(:,2)-l2(:,3))];
% for i = 1:3
%     ew(:,i) = (e0(:,i)+e1(:,i)+e2(:,i))/3;
% end

s=(de0+de1+de2)/2;
%Af - face area vector
Af=sqrt(abs(s.*(s-de0).*(s-de1).*(s-de2)));%herons formula for triangle area, 

Nverts = size(FV.vertices,1);
Nfaces = size(FV.faces,1);

%calculate weights
Acorner=zeros(Nfaces,3);
Avertex=zeros(Nverts,1);

% Calculate Vertice Normals
VertexNormals=zeros([Nverts 3]);
up=zeros([Nverts 3]);
vp=zeros([Nverts 3]);

    wfv1=Af./(de1.^2.*de2.^2);
    wfv2=Af./(de0.^2.*de2.^2);
    wfv3=Af./(de1.^2.*de0.^2);

for i=1:Nfaces

%     Calculate weights according to N.Max [1999]
    
%     wfv1=Af(i)/(de1(i)^2*de2(i)^2);
%     wfv2=Af(i)/(de0(i)^2*de2(i)^2);
%     wfv3=Af(i)/(de1(i)^2*de0(i)^2);
    
    VertexNormals(FV.faces(i,1),:)=VertexNormals(FV.faces(i,1),:)+wfv1(i)*N(i,:);
    VertexNormals(FV.faces(i,2),:)=VertexNormals(FV.faces(i,2),:)+wfv2(i)*N(i,:);
    VertexNormals(FV.faces(i,3),:)=VertexNormals(FV.faces(i,3),:)+wfv3(i)*N(i,:);
    %Calculate areas for weights according to Meyer et al. [2002]
    %check if the tringle is obtuse, right or acute
    
    if ew(i,1)<=0
        Acorner(i,2)=-0.25*l2(i,3)*Af(i)/(e0(i,:)*e2(i,:)');
        Acorner(i,3)=-0.25*l2(i,2)*Af(i)/(e0(i,:)*e1(i,:)');
        Acorner(i,1)=Af(i)-Acorner(i,2)-Acorner(i,3);
    elseif ew(i,2)<=0
        Acorner(i,3)=-0.25*l2(i,1)*Af(i)/(e1(i,:)*e0(i,:)');
        Acorner(i,1)=-0.25*l2(i,3)*Af(i)/(e1(i,:)*e2(i,:)');
        Acorner(i,2)=Af(i)-Acorner(i,1)-Acorner(i,3);
    elseif ew(i,3)<=0
        Acorner(i,1)=-0.25*l2(i,2)*Af(i)/(e2(i,:)*e1(i,:)');
        Acorner(i,2)=-0.25*l2(i,1)*Af(i)/(e2(i,:)*e0(i,:)');
        Acorner(i,3)=Af(i)-Acorner(i,1)-Acorner(i,2);
    else
        ewscale=0.5*Af(i)/(ew(i,1)+ew(i,2)+ew(i,3));
        Acorner(i,1)=ewscale*(ew(i,2)+ew(i,3));
        Acorner(i,2)=ewscale*(ew(i,1)+ew(i,3));
        Acorner(i,3)=ewscale*(ew(i,2)+ew(i,1));
    end
    Avertex(FV.faces(i,1))=Avertex(FV.faces(i,1))+Acorner(i,1);
    Avertex(FV.faces(i,2))=Avertex(FV.faces(i,2))+Acorner(i,2);
    Avertex(FV.faces(i,3))=Avertex(FV.faces(i,3))+Acorner(i,3);
    
%     Calculate initial coordinate system
    up(FV.faces(i,1),:)=e2_norm(i,:);
    up(FV.faces(i,2),:)=e0_norm(i,:);
    up(FV.faces(i,3),:)=e1_norm(i,:);
end


% for i=1:1:Nverts
%     find(FV.faces==1)
%     VertexNormals(i,1:3) = N()
%     
%     
%     up(FV.faces(i,1),:)=e2_norm(i,:);
%     up(FV.faces(i,2),:)=e0_norm(i,:);
%     up(FV.faces(i,3),:)=e1_norm(i,:);
% end

VertexNormals=normr(VertexNormals);

%Calculate initial vertex coordinate system
up = normr(cross(up,VertexNormals));
vp = cross(VertexNormals,up);

% for i=1:Nverts
%     up(i,:)=cross(up(i,:),VertexNormals(i,:));
%     up(i,:)=up(i,:)/norm(up(i,:));
%     boh(i,1) = norm(up(i,:));
%     vp(i,:)=cross(VertexNormals(i,:),up(i,:));
% end

% disp('Finished Calculating vertex normals');
end
