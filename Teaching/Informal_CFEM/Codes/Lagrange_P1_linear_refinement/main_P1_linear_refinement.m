clearvars; close all;

%----------------------
% Two point elliptic boundary value problem
%
% ODE:  -d^2y/dx^2 = f        in (0,1)
%  BC:    y(0) = 0 = y(1)
%-----------------------


%%  Domain 
x0 = 0; xfin  = 1;   % initial and final spatial nodes
N = 4; M = (2^N) +1;
% The grid consists of 2^N + 1 points. This will ensure
% that the mesh discretisation parameter 
% h(N) = mesh(i+1) - mesh(i) = 2^-N
mesh  = linspace(0,1,M); % 2^N grid points

h = mesh(2) - mesh(1);  % discretisation parameter

%% Element structure
% # (elements) = 2^N  = M-1
% # (nodes) = 2^N + 1 = M

%% Initialization 
A = sparse(M,M); % Coefficient matrix
b = sparse(M,1); % Load vector

%% Assembly 
for nelem = 1:M-1
    A(nelem:nelem+1,nelem:nelem+1) = A(nelem:nelem+1,nelem:nelem+1) +...
                                       (1/h)*[1 -1;-1 1];
    % A(i1:i2,j1:j2) will extract the block in A from rows i1 to i2 and 
    % columns j1 to j2. The last expression extracts the block of A
    % corresponding to element # (nlem) (i.e, rows nlem to nlem+1 and
    % column nlem to nlem + 1) and adds the *local stiffness matrix*
    % [1/h -1/h;-1/h 1/h] to it. 
    xmid = mesh(nelem) + h/2;   % midpoint of element # (nlem)
    b(nelem:nelem+1,1) = b(nelem:nelem+1,1) +...
                             (h/2)*load_new(xmid)*[1;1];
end

%% Linear solving 
bd_var = [1;M];  % Boundary variables with known values
all_var = [1:M]'; % All variables 

free_var = setdiff(all_var,bd_var);  % Unknown variables

soln = zeros(M,1);     
% soln in the M+1*1 column vector to store the nodal 
% values of y_h. In the next step, we are assigning 
% boundary values y(0) = 0 and y(1) = 0 to the 
% boundary nodes.
soln(1) = 0; soln(M) = 0;
%soln(bd_var) = 0;

soln(free_var) = A(free_var,free_var)\b(free_var);

%% Plotting 
% Exact solution
figure(1)
plot(mesh,uex(mesh),'blue','LineWidth',1.5);
title('Exact solution','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y_{ex}$','Interpreter','latex');

% Approximate solution
figure(2)
plot(mesh,soln,'red','LineWidth',1.5);
title('Approximate solution','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y_{h}$','Interpreter','latex');







