function [f,dfdx] = vanderpol_jac(t,x,nu)
% Van der Pol system with Jacobian

f = [x(2) ; 
     nu*(1-x(1)^2)*x(2) - x(1)  ];

dfdx = [ 0 1;
         -nu*2*x(1)*x(2)-1 nu*(1-x(1)^2)];

