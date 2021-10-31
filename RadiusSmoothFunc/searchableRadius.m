function [SearchPosInd] = searchableRadius(sphCent, V, SearchRadius, Nsmooth)
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
%% [SearchPosInd] = searchRadius(sphCent, V, SearchRadius)
% Find sorted Indexes on the {V}  indexing reference that are within or on the surface the sphere 
% built on the point {sphCent} with radius [SearchRadius]. The indexes are sorted from the farthest 
% to the closest to the {sphCent}. 
sphCent = repmat(sphCent,length(V),1);
% V = cell2mat(V);
% SearchRadius = cell2mat(SearchRadius);
SearchPointsInd = find(sum(V <= sphCent+SearchRadius & V >= sphCent-SearchRadius,2)==3); % selected searchable indexes on original V reference
SearchPointsXYZ = V(sum(V <= sphCent+SearchRadius & V >= sphCent-SearchRadius,2)==3,:);  % selected searchable points
if size(SearchPointsXYZ,1) <= Nsmooth
    SearchPointsInd = find(sum(V <= sphCent+SearchRadius*10 & V >= sphCent-SearchRadius*10,2)==3); % selected searchable indexes on original V reference
    SearchPointsXYZ = V(sum(V <= sphCent+SearchRadius*10 & V >= sphCent-SearchRadius*10,2)==3,:);  % selected searchable points
end
    
Dxyz = sphCent(1:size(SearchPointsXYZ,1),:)-SearchPointsXYZ; % distances x,y,z
SearchPointsInd(:,2) = sqrt(Dxyz(:,1).^2+Dxyz(:,2).^2+Dxyz(:,3).^2); % distance (~radius)

% sorting the indexes from the highest to the lowest (for using an exponential moving average purpose)
SearchPointsInd = sortrows(SearchPointsInd,-2);

% sorted Indexes on the V reference that are within the sphere built on the point [sphCent] with radius [SearchRadius]
SearchPosInd = single(SearchPointsInd(SearchPointsInd(:,2)<=SearchRadius*10,1));
Nind = min([Nsmooth+3 length(SearchPosInd)-1]); 
SearchPosInd = SearchPosInd(end-Nind:end,1);

end