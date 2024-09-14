function out  = error_H1(elem,coord,soln,el2)
% This code computes the H1 error between the exact and the approximate
% solution. Please see the comments for the L2 error computation. 

% error_H1 = (\int (du_approx/dx - du_exct/dx)^2dx + error_L2^2)^{1/2}

current_error = 0; 
maxelems = size(elem,1);

for j = 1:maxelems
   current_coord = coord(elem(j,:),:);
   xmid  = sum(current_coord)*0.5;
   h = abs(current_coord(1) - current_coord(2));
   diff_val = (soln(elem(j,2)) - soln(elem(j,1)))/h;
   local_error  = h*(diff_val - duex(xmid))^2;
   current_error = current_error + local_error;
end

  out = sqrt(current_error + el2);
end
