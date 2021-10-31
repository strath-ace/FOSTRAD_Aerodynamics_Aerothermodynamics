function [ A ] = Area( vert_Tri, points )
%UNTITLED5 Summary of this function goes here
%   A = Area( faces, verts )


for i = 1 : length(vert_Tri)
    
V1 = vert_Tri(i,1);
V2 = vert_Tri(i,2);
V3 = vert_Tri(i,3);

% A = tri_area(P1, P2, P3)
%
% DESC:
% calculates the triangle area given the triangle vertices (using Heron's
% formula)
%
% AUTHOR
% Marco Zuliani - zuliani@ece.ucsb.edu
%
% VERSION:
% 1.0
%
% INPUT:
% P1, P2, P3 = triangle vertices
%
% OUTPUT:
% A          = triangle area

P1 = points(V1,:);
P2 = points(V2,:);
P3 = points(V3,:);

u1 = P1 - P2;
u2 = P1 - P3;
u3 = P3 - P2;

a = norm(u1);
b = norm(u2);
c = norm(u3);

% Stabilized Heron formula
% see: http://http.cs.berkeley.edu/%7Ewkahan/Triangle.pdf
%
% s = semiperimeter
% A = sqrt(s * (s-a) * (s-b) * (s-c))

% sort the elements
v = sort([a b c]);
a = v(3);
b = v(2);
c = v(1);

temp = b + c;
v1 = a + temp;
temp = a - b;
v2 = c - temp;
v3 = c + temp;
temp = b - c;
v4 = a + temp;

A(i,1) = 0.25 * sqrt(abs(v1*v2*v3*v4));



end

end

