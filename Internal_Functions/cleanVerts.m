%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOSTRAD: The Free Open Source Tool for Re-entry of Asteroids and Debris %%
%%%%%%%%  Copyright (C) 2018 University of Strathclyde and Authors %%%%%%%%%%
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

function [varargout] = cleanVerts(F,V,N)
    % [F1, V1, N1] = cleanVerts(F,V,N)
    % removes unused vertexes

    
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