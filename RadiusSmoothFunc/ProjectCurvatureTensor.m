function [new_ku,new_kuv,new_kv]=ProjectCurvatureTensor(uf,vf,nf,old_ku,old_kuv,old_kv,up,vp)
%{ Summary: ProjectCurvatureTensor performs a projection
%of the tensor variables to the vertexcoordinate system
%INPUT:
%uf,vf - face coordinate system
%old_ku,old_kuv,old_kv - face curvature tensor variables
%up,vp - vertex cordinate system
%OUTPUT:
%new_ku,new_kuv,new_kv - vertex curvature tensor variabels
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
[r_new_u,r_new_v]=RotateCoordinateSystem(up,vp,nf);
OldTensor=[old_ku old_kuv; old_kuv old_kv];
u1=r_new_u*uf;
v1=r_new_u*vf;
u2=r_new_v*uf;
v2=r_new_v*vf;
new_ku=[u1 v1]*OldTensor*[u1;v1];
new_kuv=[u1 v1]*OldTensor*[u2;v2];
new_kv=[u2 v2]*OldTensor*[u2;v2];
%{
new_ku=old_ku*u1*u1+2*old_kuv*(u1*v1)+old_kv*v1*v1;
new_kuv=old_ku*u1*u2+old_kuv*(u1*v2+u2*v1)+old_kv*v1*v2;
new_kv=old_ku*u2*u2+2*old_kuv*(u2*v2)+old_kv*v2*v2;
%}
end
