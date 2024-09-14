clearvars; close all; clc;
format long;

%---------------------------------------------%
% Two point elliptic boundary value problem   %
% with red-refinement using sparse assembly   %
%                                             %
% ODE:  -d^2y/dx^2 = f        in (0,1)        %
%  BC:    y(0) = a  y(1) = b                  %
%          %
%---------------------------------------------%


%%  Call initial geometry
% The function initial_geometry(n) will load two initial geometries for the
% case n = 1 and n = 2. The outputs are elements and coordinates. 
% n = 1:  elem = [1 2;2 3] coord = [0;0.5;1]
% n = 2:  elem = [1 2;2 3;3 4] coord = [0;1/3;2/3;1]
[elem,coord] = initial_geometry(1);

% Define the boundary nodes. The boundary node indices are not changed by
% the refinement procedure. Therefore, the boundary nodes are always the
% first and the last members of the coord vector.
bd_nodes = [1,size(coord,1)]; 

% Maximum number of refinements
maxref = 10;

% storage matrices
hval = zeros(maxref,1);    % Discretisation factor
error_val_L2 = hval;       % Error vector for L2
error_val_H1 = hval;       % Error vector for H1
L2_rate = hval;            % Rate of convergence WRT L2 error
H1_rate = hval;            % Rate of convergence WRT H1 error


%% Refinement loop
for ref = 1:maxref
    % refinement procedure: On each call of redrefine(.,.), the current
    % mesh refines with the subsimplices formed by conjoining endpoints and
    % midpoint of an interval are added as the new elements. 
    %
    %   *--------------*   : old element
    %   *------*-------*   : two new elements in the refined mesh
[elem,coord] = redrefine(elem,coord);    

% Data from element-coord structures
M = size(coord,1);       % # (nodes/coordinates) = M
maxelems = size(elem,1); % # (elements)
h = abs(coord(elem(1,1)) - coord(elem(1,2)));  % discretisation factor
hval(ref) = h;    % store h the storage matrix 


%% Sparse assembly

% Sparse allocation of memory to the coefficient and load vector. Sparsity
% saves memory.

tic
row_address = reshape(repmat(reshape(elem',size(elem,1)*size(elem,2),1),...
                 1,2)',size(elem,1)*size(elem,2)^2,1);
column_address = reshape(repmat(elem,1,size(elem,2))',1,size(elem,1)*size(elem,2)^2);
element_entries = (1/h)*repmat([1 -1 -1 1],1,size(elem,1));
A = sparse(row_address,column_address,element_entries,M,M);

f_mid  = load_new(0.5*sum(coord(elem),2));
cm_f = reshape(elem',size(elem,1)*size(elem,2),1);
v_f = reshape(repmat(f_mid,1,2)',size(elem,1)*size(elem,2),1);
b = (h/2)*sparse(cm_f,ones(size(elem,1)*size(elem,2),1),v_f,M,1);
toc

%% Linear solving
all_var = [1:M]'; % All variables 

% Since boundary data is available, that needs to be removed from all
% variables. The free variables are then given by:
free_var = setdiff(all_var,bd_nodes);  

soln = zeros(M,1);     
% soln in the M*1 column vector to store the nodal  values of y_h. In the
% next step, we are assigning boundary values y(0) = 0 and y(1) = 0 to the 
% boundary nodes.
% soln(bd_nodes) = 0;
soln(bd_nodes) = uex(coord(bd_nodes,:));
bnew = b - A*soln;

% The solving opration is presented in the next step. 
% soln(at free variables) = A(rows and columns corresponding to free
%                                   variables)^-1 b(entries corresponding 
%                                                   to the free variables)
soln(free_var) = A(free_var,free_var)\bnew(free_var);

% This step calculates the L2 and H1 errors between the approximate 
% solution and the exact solution.
error_val_L2(ref) = error_L2(elem,coord,soln);
error_val_H1(ref) = error_H1(elem,coord,soln,error_val_L2(ref));
end


% Computation of asymptotic convergence rates
for jref = 1:maxref-1
    L2_rate(jref+1) = log(error_val_L2(jref)/error_val_L2(jref+1))/...
                           log(hval(jref)/hval(jref+1));
    H1_rate(jref+1) = log(error_val_H1(jref)/error_val_H1(jref+1))/...
                           log(hval(jref)/hval(jref+1));
end

table(hval,error_val_L2,L2_rate,error_val_H1,H1_rate)

%% Graphical represenation of the extact and approximate solution
% Exact solution
figure(1)
for i4plot = 1:maxelems
    cur_coord = coord(elem(i4plot,:),:);
    uval = uex(cur_coord);
    plot(cur_coord,uval,'blue','LineWidth',1.5);
    hold on;
end
title('Exact solution','Interpreter','latex','FontSize',14);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$y_{ex}$','Interpreter','latex','FontSize',14);
hold off;

% Approximate solution
figure(2)
for i4plot = 1:maxelems
    cur_coord = coord(elem(i4plot,:),:);
    uval = soln(elem(i4plot,:),:);
    plot(cur_coord,uval,'blue','LineWidth',1.5);
    hold on;
end
title('Approximate solution','Interpreter','latex','FontSize',14);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$y_{h}$','Interpreter','latex','FontSize',14);







