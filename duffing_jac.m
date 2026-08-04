function [f,dfdx] = duffing_jac(t,x,alpha,beta,gamma,delta,omega)
% Duffing system with Jacobian

f = [x(2) ; 
     -delta*x(2)-alpha*x(1)-beta*x(1)^3+gamma*cos(omega*t) ];

dfdx = [ 0 1;
         -alpha-3*beta*x(1)^2 -delta];

