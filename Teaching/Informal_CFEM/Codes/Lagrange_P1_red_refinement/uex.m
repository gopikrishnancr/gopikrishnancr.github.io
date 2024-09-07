function out = uex(x)
out  = x.*(1-x);
% .* is the elementwise multiplication operator. This
% overloads the function defintion to accept vector
% inputs also.
