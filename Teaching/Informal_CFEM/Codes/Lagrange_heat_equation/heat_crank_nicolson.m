clearvars; close all; clc;
format long;

%---------------------------------------------%
% Two point elliptic boundary value problem   %
% with red-refinement using sparse assembly   %
%                                             %
% ODE:  du/dt - d^2u/dx^2 = f in (a,b)        %
%  BC:    u(t,a) = u_a(t)  u(t,b) = u_b(t)    %
%  IC:    u(0,x) = u_0(x)                     %
%---------------------------------------------%


%%  Call initial geometry
% The function initial_geometry(n) will load two initial geometries for the
% case n = 1 and n = 2. The outputs are elements and coordinates. 
% n = 1:  elem = [1 2;2 3] coord = [0;0.5;1]
% n = 2:  elem = [1 2;2 3;3 4] coord = [0;1/3;2/3;1]
[elem,coord] = initial_geometry(1);

% Define the boundary nodes. The boundary node indices are not changed by
% the refinement procedure. Therefore, the boundary nodes are always the
% first and the third members of the coord vector.
bd_nodes = [1,3]; 

% Maximum number of refinements
maxref = 4;

% storage matrices
hval = zeros(maxref,1);    % Discretisation factor
error_val_L2 = hval;       % Error vector for L2
error_val_H1 = hval;       % Error vector for H1
L2_rate = hval;            % Rate of convergence WRT L2 error
H1_rate = hval;            % Rate of convergence WRT H1 error

% Crank-Nicolson parameter 
theta = 0;

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
ncoord = size(coord,1);       % # (nodes/coordinates) = M
maxelems = size(elem,1); % # (elements)
h = abs(coord(elem(1,1)) - coord(elem(1,2)));  % discretisation factor
hval(ref) = h;    % store h the storage matrix 

k = 0.1*h^2;         % CFL condition
tmesh  = [0:k:1]';     % time-mesh

%% Time independent construction 
%       Sparse assembly
row_address = reshape(repmat(reshape(elem',size(elem,1)*size(elem,2),1),...
                 1,2)',size(elem,1)*size(elem,2)^2,1);
column_address = reshape(repmat(elem,1,size(elem,2))',1,size(elem,1)*size(elem,2)^2);
element_entries_stiffness = (1/h)*repmat([1 -1 -1 1],1,size(elem,1));
S = sparse(row_address,column_address,element_entries_stiffness,ncoord,ncoord);
element_entries_mass = h*repmat([1/3 1/6 1/6 1/3],1,size(elem,1));
M  =  sparse(row_address,column_address,element_entries_mass,ncoord,ncoord);

%      Boundary data
all_var = [1:ncoord]'; % All variables 
free_var = setdiff(all_var,bd_nodes);

%      Storage matrix for approximate solution
uval  = sparse(ncoord,length(tmesh));

%      Crank-Nicolson matrices
R = M + k*theta*S;   L  = M - k*(1-theta)*S;

%      Address for load vector
cm_f = reshape(elem',size(elem,1)*size(elem,2),1);

%      Initial value 
uval(:,1) = heat_exact(0,coord);

 for tn = 1:length(tmesh)-1
        % Assembly of the time-dependent load term
        CN_time = theta*tmesh(tn+1) + (1-theta)*tmesh(tn);
        f_mid  = heat_load(CN_time,0.5*sum(coord(elem),2));
        v_f = reshape(repmat(f_mid,1,2)',size(elem,1)*size(elem,2),1);
        b = (h/2)*sparse(cm_f,ones(size(elem,1)*size(elem,2),1),v_f,ncoord,1);

        % Impose the boundary data
        uval(bd_nodes,tn+1) = heat_exact(tmesh(tn+1),coord(bd_nodes));
        bnew = (L*uval(:,tn) + k*b) - R*uval(:,tn+1);
        uval(free_var,tn+1) = R(free_var,free_var)\bnew(free_var);
 end

error_val_L2(ref) = error_L2(elem,coord,uval(:,end));
error_val_H1(ref) = error_H1(elem,coord,uval(:,end),error_val_L2(ref));
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
    uval_current = heat_exact(1,cur_coord);
    plot(cur_coord,uval_current,'blue','LineWidth',1.5);
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
    uval_cur = uval(elem(i4plot,:),end);
    plot(cur_coord,uval_cur,'blue','LineWidth',1.5);
    hold on;
end
title('Approximate solution','Interpreter','latex','FontSize',14);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$y_{h}$','Interpreter','latex','FontSize',14);







