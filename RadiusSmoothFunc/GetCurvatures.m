function [CurvUxVy,Cmagnitude, Avertex, Acorner, normalVerts , avEdge] = GetCurvatures( F, V, N, toggleDerivatives);
% function [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude] = GetCurvatures( FV,toggleDerivatives )
% Summary 
%Author: Itzik Ben Shabat
%Last Update: July 2014
%Implemented according to "Estimating Curvatures and Their Derivatives on Triangle Meshes" by Szymon Rusinkiewicz (2004)
%and according to its C implementation trimesh2
%{
%GetCurvatures computes the curvature tensor and the principal curvtures at
%each vertex of a mesh given in a face vertex data structure
%INPUT: 
%FV -struct - Triangle mesh face vertex data structure (containing FV.face and
%FV.Vertices) 
%N - triangles mesh normals on faces
%toggleDerivatives - scalar  1 or 0 indicating wether or not to calcualte curvature derivatives
%OUTPUT:
% PrincipalCurvatures - 2XN matrix (where N is the number of vertices
%       containing the proncipal curvtures k1 and k2 at each vertex
% PrincipalDir1 - NX3 matrix containing the direction of the k1 principal
%       curvature
% PrincipalDir2 - NX3 matrix containing the direction of the k2 principal
%       curvature
% FaceCMatrix - 4XM matrix (where M is the number of faces) containing the 4
%       coefficients of the curvature tensr of each face
% VertexCMatrix- 4XN matrix (where M is the number of faces) containing the 4
%       coefficients of the curvature tensr of each tensor
% Cmagnitude - NX1 matrix containing the square sum of the curvature tensor coefficients at each
%       vertex (invariant scalar indicating the change of curvature)
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
%%
if ~exist('N','var')
    error('The new version requires the input of normals on faces')
end
tstart = tic;
FV.faces = F;
FV.vertices = V;
FaceCMatrix=NaN;
VertexCMatrix=NaN;
Cmagnitude=NaN;

% [normalFaces]=CalcFaceNormals(FV);
[normalVerts,Avertex,Acorner,up,vp,avEdge]=CalcVertexNormals(FV,N);
[FaceSFM,VertexSFM,wfp]=CalcCurvature(FV,normalVerts,N,Avertex,Acorner,up,vp);
[CurvUxVy,PrincipalDir1,PrincipalDir2]=getPrincipalCurvatures(FV,VertexSFM,up,vp);
if toggleDerivatives==1
    [FaceCMatrix,VertexCMatrix,Cmagnitude]=CalcCurvatureDerivative(FV,N,CurvUxVy,up,vp,wfp);
end
disp(['Surface Curvature Defined in: ',num2str(toc(tstart)),'s'])
end

