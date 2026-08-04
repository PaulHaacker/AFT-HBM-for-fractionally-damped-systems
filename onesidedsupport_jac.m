function [f,dfdx] = onesidedsupport_jac(t,x,alpha,beta,gamma,delta,omega)
% one-DOF oscillator with one-sided support

f = [x(2) ; 
     -delta*x(2)-alpha*x(1)-beta*max(x(1),0)+gamma*cos(omega*t) ];

if x(1)>0
 dsigmadx = 1;
else
 dsigmadx = 0;
end

dfdx = [ 0 1;
         -alpha-beta*dsigmadx -delta];


