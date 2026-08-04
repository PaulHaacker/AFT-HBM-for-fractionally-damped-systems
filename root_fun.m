function [f,dfdx] = root_fun(x)
% two-dimensional root-finding problem
f = [x(1) + x(2)^2 + 1;
     x(2) + x(2)^3 + x(1)^2-3];

dfdx = [1 2*x(2);
        2*x(1) 1+3*x(2)^2];