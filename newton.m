function [x,f,dfdx,conv,k] = newton(fun,x0,tol,maxiter)
% Newton's method to find zeros
%
% Input
%  fun     : function which gives f and dfdx
%  x0      : initial guess
%  tol     : tolerance
%  maxiter : maximum number of iterations
%
% Output
%  x       : solution for which f(x) = 0 within tolerance
%  f       : function value f(x)
%  dfdx    : Jacobian matrix
%  conv    : boolean to indicate convergence
%  k       : number of iterations

% Remco Leine, INM, University of Stuttgart, 2023
x = x0;
k = 0;
conv = 0;
while ~conv && k<maxiter
    [f,dfdx] = fun(x);
    if norm(f) < tol
        conv = 1;
    else    
        k = k+1;
        x = x - dfdx\f;
    end
end

