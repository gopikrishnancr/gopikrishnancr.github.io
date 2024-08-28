% This is a basic code tha implements the Lagrange P1 FEM for the problem
% -d^y + y = f with initial and boundary conditions 
% y(0) = 0 and y(1) = 0.
%-------------------------------------------------------------------------

clear all;  % clear all the variables 
% alternatively % clearvars

close all;  % close all previous windows

format LONGE;  % to display numbers in standad scientific way

% Domain or interval 
xini = 0; xfin  = 1;  N = 5;

xmesh  = linspace(xini,xfin,N)';    % ' transpose opr
h = xmesh(2) - xmesh(1);            % mesh parameter

stimat = (1/h)*[1 -1;-1 1]; % stiffness matrix
massmat  = h*[1/3 1/6;1/6 1/3]; % mass matrix
A = sparse(N,N); b = zeros(N,1);


for j = 1:N-1
    A(j:j+1,j:j+1) = A(j:j+1,j:j+1) + stimat + massmat;
    xmid = (xmesh(j) + xmesh(j+1))/2;
    b(j:j+1,1) =  b(j:j+1,1) + (h/2)*[f(xmid);f(xmid)];
end



uh = zeros(N,1);
uh(2:N-1) = A(2:N-1,2:N-1)\b(2:N-1);  % A^-1*b = A\b

figure(1)
plot(xmesh,uh);
hold on;
plot(xmesh,uexact(xmesh),'--','Color','red');







