% filename: redrefine.m

function  [elem_n,coord_n] = redrefine(elem,coord)

% Input:  elem - nelem * 2 matrix, nlem = #(elements)
%         coord - ncoord * 1 vector, ncoord = #(coordinates)

% The new coordinate vector is initialised as the old coordinate vector.
% As and when the mesh is refined, we shall add the new coordinates as new
% rows to coord. 
coord_n = coord;      
% elem_n is the empty array to store the new elements. 
elem_n = [];

nelem  = size(elem,1);  % # of elements
ncoord = size(coord,1); % # of coordinates

% Red-refinement procedure.
%
%

for i = 1:nelem
    % Compute the midpoint of the current element. 
    xnew = sum(coord(elem(i,:),:))*0.5;
    % This midpoint is appended as a new coordinate to the coord_n. Recall
    % that coord_n already contains the old coordinates coord (see line 11) 
    coord_n = [coord_n;xnew];
    % The new elements are 
    %      [x_i x_mid] and [xmid x_{i+1}]
    % where x_i and x_{i+1} are the left and right end points of the
    % current element. The index of x_i and x_{i+1} are elem(i,1) and
    % elem(i,2), and that of x_mid is ncoord+i (the new coordinate is after
    % i columns from the last column (ncoord^th) of old coord). 
    new_elems = [elem(i,1) ncoord + i;
                 ncoord + i elem(i,2)];
    elem_n = [elem_n;new_elems];
end
