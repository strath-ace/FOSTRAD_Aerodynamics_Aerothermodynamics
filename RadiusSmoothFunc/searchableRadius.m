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

function [SearchPosInd] = searchableRadius(sphCent, V, SearchRadius, Nsmooth)
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