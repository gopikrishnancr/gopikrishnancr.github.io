function out = uex(x)
out  = x.^2 + sin(x) + 2;
% .* is the elementwise multiplication operator. This
% overloads the function defintion to accept vector
% inputs also.
