function [varargout] = cleanVerts(F,V,N)
    % [F1, V1, N1] = cleanVerts(F,V,N)
    % removes unused vertexes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    
    F1 = reshape(F,[],1);
    if exist('N','var')
        N1 = reshape(N,[],1);
    end
    
    Vtemp = [V(F,:) F1];
    [V1, ~, indV1] = unique(round(Vtemp(:,1:3),8),'rows','stable');
    % V1 = V1(:,1:end-1);
    F1 = reshape(indV1,[],size(F,2));
    if exist('N','var')
        N1 = reshape(N(indV1),[],3);
    end
    
    if nargout == 3 && nargin == 3
        varargout = {F1 V1 N1};
    else
        varargout = {F1 V1};
    end
end