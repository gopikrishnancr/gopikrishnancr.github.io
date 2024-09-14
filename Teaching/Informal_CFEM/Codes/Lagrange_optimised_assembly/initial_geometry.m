% Initial Geometry 
% filename: initial_geometry.m

function [Elem,Coord] = initial_geometry(n)

switch n
    case 1
        Elem = [1 2;2 3];
        Coord = [0; 0.5; 1];
    case 2
        Elem = [1 2;2 3;3 4];
        Coord = [0;1/3;2/3;1];
end
