clearvars; close all;

%----------------------
% Two point elliptic boundary value problem
%
% ODE:  dy/dt - d^2y/dx^2 = f  in (0,1)
%  BC:    y(t,0) = a(t)   y(t,1) = b(t)
%  IC:    y(0,x) = y0(x)
%-----------------------


%%  Domain 
x0 = 0; xfin  = 1;   % initial and final spatial nodes
N = 5; M = (2^N) +1;
% The grid consists of 2^N + 1 points. This will ensure
% that the mesh discretisation parameter 
% h(N) = mesh(i+1) - mesh(i) = 2^-N
mesh  = linspace(0,1,M)'; % 2^N grid points

h = mesh(2) - mesh(1);  % discretisation parameter

% CFL condition: k/h^2 < 1 - Ensure stability of the scheme
% k = 0.5*h^2

k = h^2;    time_mesh  = [0:k:1]';

%% Element structure
% # (elements) = 2^N  = M-1
% # (nodes) = 2^N + 1 = M

%% Initialization 
Mass_mat = sparse(M,M); % Mass matrix
Stiff_mat = sparse(M,M);

yval  = sparse(M,length(time_mesh));
% yval(:,n) - y(t_n,.) - approximate solution at the nth time step

yval(:,1) = uex(0,mesh);   % Initial condition



%% Assembly 
for nelem = 1:M-1
    Stiff_mat(nelem:nelem+1,nelem:nelem+1) = Stiff_mat(nelem:nelem+1,nelem:nelem+1) +...
                                       (1/h)*[1 -1;-1 1];
    Mass_mat(nelem:nelem+1,nelem:nelem+1) = Mass_mat(nelem:nelem+1,nelem:nelem+1) +...
                                       (h)*[1/3 1/6;1/6 1/3];
end
 
% Time mesh
  % midpoint of element # (nlem)
bd_var = [1;M];  % Boundary variables with known values
all_var = [1:M]'; % All variables 
free_var = setdiff(all_var,bd_var);  % Unknown variables

for tn = 1:length(time_mesh)-1
    b = sparse(M,1); % Load vector
    for nelem = 1:M-1
    xmid = mesh(nelem) + h/2; 
    b(nelem:nelem+1,1) = b(nelem:nelem+1,1) +...
                             (h/2)*load_new(time_mesh(tn),xmid)*[1;1];
    end
    F = (Mass_mat - k*Stiff_mat)*yval(:,tn) + k*b;
    yval(1,tn+1) = uex(time_mesh(tn+1),0); 
    yval(M,tn+1) = uex(time_mesh(tn+1),1); 

    F = F - M*yval(:,tn+1);
    yval(free_var,tn+1) = Mass_mat(free_var,free_var)\F(free_var);
end


%% Plotting 
% Exact solution
% figure(1)
% plot(mesh,uex(mesh),'blue','LineWidth',1.5);
% title('Exact solution','Interpreter','latex');
% xlabel('$x$','Interpreter','latex');
% ylabel('$y_{ex}$','Interpreter','latex');

% Approximate solution
figure(2)
for tn = 1:length(time_mesh)
plot(mesh,yval(:,tn),'red','LineWidth',1.5);
title('Approximate solution','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y_{h}$','Interpreter','latex');
axis square;
xlim([0 1]); ylim([0 1]);
pause(0.001);
end






