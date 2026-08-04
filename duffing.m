function f = duffing(t,x,alpha,beta,gamma,delta,omega)
% Duffing system 

f = [x(2) ; 
     -delta*x(2)-alpha*x(1)-beta*x(1)^3+gamma*cos(omega*t) ];


