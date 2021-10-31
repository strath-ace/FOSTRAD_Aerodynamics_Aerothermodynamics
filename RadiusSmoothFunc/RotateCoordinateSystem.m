function [r_new_u,r_new_v]=RotateCoordinateSystem(up,vp,nf)
%{Summary: RotateCoordinateSystem performs the rotation of the vectors up and vp
%to the plane defined by nf
%INPUT:
%up,vp- vectors to be rotated (vertex coordinate system)
%nf - face normal
%OUTPUT:
%r_new_u,r_new_v - rotated coordinate system
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
r_new_u=up;
r_new_v=vp;
np=cross(up,vp);
np=np/norm(np);
ndot=nf*np';
if ndot<=-1
    r_new_u=-r_new_u;
    r_new_v=-r_new_v;  
    return;
end
perp=nf-ndot*np;
dperp=(np+nf)/(1+ndot);
r_new_u=r_new_u-dperp*(perp*r_new_u');
r_new_v=r_new_v-dperp*(perp*r_new_v');
end