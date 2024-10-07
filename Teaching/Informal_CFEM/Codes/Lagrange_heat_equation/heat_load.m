function out = heat_load(t,x)
out  = exp(-t).*sin(pi*x).*(pi^2 - 1);
