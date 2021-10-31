function [PrincipalCurvatures,PrincipalDir1,PrincipalDir2]=getPrincipalCurvatures(FV,VertexSFM,up,vp);
%% Summary 
%Author: Itzik Ben Shabat
%Last Update: July 2014
%{
%getPrincipalCurvatures calculates the pincipal curvatures and their
%directions from the second fundemental tensor
%INPUT: FV -  mesh face-vertex data structure
%VertexSFM - second fundemental matrix for each vertex
%up,vp - vertex local coordinate frame
%OUTPUT:
%PrincipalCurvatures - 2X[number of vertices] matrix containing the
%principal curvature values
%PrincipalDir1,PrincipalDir2 - direction vectors of the proncipal
%components
%}
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
%
%% Code
% disp(' Calculating Principal Components');

%Calculate Principal Curvatures
PrincipalCurvatures=zeros(size(FV.vertices,1),1);
PrincipalDir1=zeros(size(FV.vertices,1),3);
PrincipalDir2=zeros(size(FV.vertices,1),3);
%{
%solve using matlab eigen value and eigenvector solver
for i=1:size(FV.vertices)
    [V,D]=eig(VertexSFM{i,1});
    % PrincipalCurvatures(:,i)=eig(VertexSFM{i,1});
    [maxk indexk]=max(abs([D(1,1),D(2,2)]));
    if indexk==1
        PrincipalCurvatures(1,i)=D(1,1);
        PrincipalCurvatures(2,i)=D(2,2);
        PrincipalDir1(i,:)=(V(1,1)*up(i,:)+V(2,1)*vp(i,:));
        PrincipalDir2(i,:)=(V(1,2)*up(i,:)+V(2,2)*vp(i,:));
    else
        PrincipalCurvatures(2,i)=D(1,1);
        PrincipalCurvatures(1,i)=D(2,2);
        PrincipalDir2(i,:)=(V(1,1)*up(i,:)+V(2,1)*vp(i,:));
        PrincipalDir1(i,:)=(V(1,2)*up(i,:)+V(2,2)*vp(i,:));
    end
end
%}
%
for i=1:size(FV.vertices)
%This is taken from trimsh2 - Szymon Rusinkiewicz implementation. It also considers
%direction
np=cross(up(i,:),vp(i,:));
[r_old_u, r_old_v]=RotateCoordinateSystem(up(i,:), vp(i,:), np);
ku=VertexSFM{i,1}(1,1);
kuv=VertexSFM{i,1}(1,2);
kv=VertexSFM{i,1}(2,2);
c = 1;
s = 0;
tt = 0;
if kuv ~= 0
    %Jacobi rotation to diagonalize
    h = 0.5 * (kv - ku) / kuv;
    if  h < 0
        tt =1 / (h - sqrt(1 + h*h)) ;
    else
        tt =1/ (h + sqrt(1 + h*h));
    end
    c = 1 / sqrt(1+ tt*tt);
    s = tt * c;
end

k1 = ku - tt * kuv;
k2 = kv + tt * kuv;

if (abs(k1) >= abs(k2))
    PrincipalDir1(i,:) = c*r_old_u - s*r_old_v;
else
    temp=k1;
    k1=k2;
    k2=temp;
    PrincipalDir1(i,:) = s*r_old_u + c*r_old_v;
end
PrincipalDir2(i,:) = cross(np , PrincipalDir1(i,:));
PrincipalCurvatures(i,1)=k1;
PrincipalCurvatures(i,2)=k2;

end

% if sum(isnan(PrincipalCurvatures))>=1
%     disp('nanValues have been set equal to the average curvature')
% end
    
% disp('Finished Calculating Principal Components');
end