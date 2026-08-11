function [f,dfdx,dfdxalpha] = onesidedsupport_jac_frac(t,x,x_alpha,mu,beta,gamma,delta,omega)
% one-DOF oscillator with one-sided support

f = [x(2) ; 
     -delta*x_alpha(1)-mu*x(1)-beta*max(x(1),0)+gamma*cos(omega*t) ];

if x(1)>0
 dsigmadx = 1;
else
 dsigmadx = 0;
end

dfdx = [ 0 1;
         -mu-beta*dsigmadx 0];

dfdxalpha = [0 0;
             -delta 0];
