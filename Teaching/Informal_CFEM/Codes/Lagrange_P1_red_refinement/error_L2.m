function out  = error_L2(elem,coord,soln)

% This subroutine compute the L2 error between the exact solution and the
% approximate solution. 

current_error = 0;           % variable to store the L2 error
maxelems = size(elem,1);

% The midpoint rule of integration is used to compute the error. 
for j = 1:maxelems
   % Compute the midpoint of the current element
   current_coord = coord(elem(j,:),:);
   xmid  = sum(current_coord)*0.5;
   h = abs(current_coord(1) - current_coord(2));

   % Evaluate the approximate solution at the midpoint. Note that since the
   % soln is linear between x_i and x_{i+1} (the coordinates corresponding
   % to the current element), the value of soln at the midpoint is simply
   % the average of soln(x_i) and soln(x_{i+1}).
   soln_mid = sum(soln(elem(j,:),:))*0.5;
    
   % calculation of local error in the current element using the midpoint
   % rule   \int_{x_i}^{x_{i+1}} fdx = (x_{i+1} - x_{i})*f(mid);
   local_error  = h*(soln_mid - uex(xmid))^2;
   current_error = current_error + local_error;
end

  out = sqrt(current_error);
end
