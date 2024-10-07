function out = uex(t,x)
out  = exp(-t)*sin(pi*x);
% .* is the elementwise multiplication operator. This
% overloads the function defintion to accept vector
% inputs also.
