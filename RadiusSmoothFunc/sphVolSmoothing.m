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

function [smoothSphProp] = sphVolSmoothing(InnerPointsInd, propOnVerts, avType, Nsmooth, MaxRefRadius, flatEdge, flatWeightFlag)

% [smoothSphProp = sphVolSmoothing(InnerPointsInd, propOnVerts, avType, Nsmooth, MaxRefRadius, flatEdge, flatWeightFlag)
% input:
% InnerPointsInd  -  {cell[Nx1]} points to be averaged (all points will be averaged) according to:
% propOnVerts     -  {cell[Nx1]} containing the properties to be smoothed e.g.: curvature, heat transfer, stresses etc..
% - avType  = char 'e' for an exponential moving average; 's' for a simple average; 
% - flatEdge= the algorithm detects the number of flatVertices, according to the ratio of flatVertices/Nsmooth within the 
%                       selected smoothing points, it assign the reference flat radius to the vertex.
%                       A vertex is considered flat when flatVertices/Nsmooth >= flatEdge.
%                   *** if flatEdge == 1, only the Vertices whose ALL nearest vertexes are "flat" will be considered flat
%                       all the remaining vertices (the ones with flatVertices/Nsmooth < 1) will be averaged.
%                   *** if flatEdge == 0, all vertices are given the flatRadiusReferenceValue, and the spherical smoothing will be
%                       DISABLED.
%                   *** if 0 < flatEdge < 1 it averages the vertices only when their flatVertices/Nsmooth < flatEdge
%                       in general, a higher flatEdge will provide a smoother result, decreasing the number of
%                       assigned flatReferenceRadius.
% - flatWeightFlag: Enables the weighting of flatVertices Reference Radius value according to the flatVertices/Nsmooth ratio on each point.
%                   *** if == 1, linear weighting of the Reference Radius 
%                   *** if == 2, sigmoid weighting of the Reference Radius (smaller edge angles, better smoothing toward the flat surfaces)
propOnVerts = propOnVerts{1}(InnerPointsInd{1});
MaxRefRadius = MaxRefRadius*0.99;

Npoints =  length(propOnVerts); 
if Npoints < Nsmooth+1
    Nsmooth = Npoints-2; 
    
end

propOnVerts = propOnVerts(end-Nsmooth:end);
NpointsNonInf =  sum(propOnVerts(end-Nsmooth:end)<MaxRefRadius);
NpointsInf = Nsmooth+1-NpointsNonInf;
    
flatRatio = (NpointsInf)/Nsmooth;
propBaseline = mean(propOnVerts(propOnVerts < MaxRefRadius));

if flatRatio >= flatEdge
    flatFlag = 1;
else 
    flatFlag = 0;
    if flatWeightFlag ==1 % linear weight
        propOnVerts(propOnVerts>=MaxRefRadius) = propBaseline+(MaxRefRadius-propBaseline)*flatRatio;
    elseif flatWeightFlag ==2 % sigmoid weight
        propOnVerts(propOnVerts>=MaxRefRadius) = propBaseline+(MaxRefRadius-propBaseline)*((1+erf(flatRatio*4-2))/2);
%         propOnVerts(propOnVerts>=MaxRefRadius) = MaxRefRadius*((1+erf(flatRatio*4-2))/2);
    end
end

switch flatFlag
    case 1
        smoothSphProp = MaxRefRadius;
    otherwise
        switch avType
            case 'e'
%                 smoothSphProp = tsmovavg(propOnVerts(propOnVerts<MaxRefRadius),avType,Nsmooth-NpointsInf,1);
                smoothSphProp = tsmovavg(propOnVerts(propOnVerts<MaxRefRadius),avType,Nsmooth-NpointsInf,1);
                smoothSphProp = smoothSphProp(end);
            case 's'
                smoothSphProp = mean(propOnVerts);
        end        
end

end
