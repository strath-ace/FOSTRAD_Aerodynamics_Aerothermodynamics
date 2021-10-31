function varargout=import_stl_fast(filename,mode)
% function to import ASCII STL files into matlab.  This update uses the
% textscan function instead of reading the file line by line, which can
% significantly reduce import times for very large STL files.
%--------------------------------------------------------------
% 
% inputs:       filename:   string 
%               mode:       integer 1 or 2, see outputs:
%
% outputs:      mode 1: [points,triangles,tri norms]
%               mode 2: [vertices, tri norms]
%
% author: Eric Trautmann, 3/31/2011 
%         etrautmann@gmail.com
%         based on code by Luigi Giaccari. 
%         future revisions will include support for binary files
%-----------------------------------------------------------------
% Ex: 
% filename = 'file.stl'     
% [p,t,tnorm] = import_stl_fast(filename,1)
%% License
% Copyright (c) 2011, Eric Trautmann
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
% * Neither the name of Physical Sciences Inc. nor the names of its
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
%%

if nargin<2
    mode=1;%default value
end

if ~(mode==1 || mode==2)
    error('invalid mode')
end

if nargout<3 && mode==1
    error('invalid input number /mode setting')
end
if nargout>2 && mode==2
    error('invalid input number /mode setting')
end

%open file
fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
if fid == -1
    error('File could not be opened, check name or path.')
end

%=====================================================

% fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.

fmt = '%*s %*s %f64 %f64 %f64 \n %*s %*s \n %*s %f64 %f64 %f64 \n %*s %f64 %f64 %f64 \n %*s %f64 %f64 %f64 \n %*s \n %*s \n';
C=textscan(fid, fmt, 'HeaderLines', 1);             % check newline format
fclose(fid);

%extract normal vectors and vertices
tnorm = cell2mat(C(1:3));
tnorm = tnorm(1:end,:); %strip off junk from last line

v1 = cell2mat(C(4:6));
v2 = cell2mat(C(7:9));
v3 = cell2mat(C(10:12));
if (C{4}(end) == NaN)
    v1 = v1(1:end-1,:); %strip off junk from last line
    v2 = v2(1:end-1,:); %strip off junk from last line
    v3 = v3(1:end-1,:); %strip off junk from last line
end

v_temp = [v1 v2 v3]';
v = zeros(3,numel(v_temp)/3);

v(:) = v_temp(:);
v = v';

    varargout = cell(1,nargout);
    switch mode
        case 1
            [p,t]=fv2pt(v,length(v)/3);%gets points and triangles

            varargout{1} = p;
            varargout{2} = t;
            varargout{3} = tnorm;
        case 2
            varargout{1} = v;
            varargout{2} = tnorm;
    end
end

%
function [p,t]=fv2pt(v,fnum)

%gets points and triangle indexes given vertex and facet number
c=size(v,1);

%triangles with vertex id data
t=zeros(3,fnum);
t(:)=1:c;

%now we have to keep unique points fro vertex
[p,i,j]=unique(v,'rows'); %now v=p(j) p(i)=v;
t(:)=j(t(:));
t=t';

end