function [f,dfdx,dfdF] = truss(x,F,a,l0,k,m,c)
% The truss system for Task 3.3
q = x(1);
qdot = x(2);

f = [qdot;
     -k/m*q*(1-l0/sqrt(a^2+q^2))-c/m*qdot+F/m];

dfdx = [ 0 1;
         -k/m*(1-l0/sqrt(a^2+q^2)+l0*q^2/(a^2+q^2)^(3/2)) -c/m];

dfdF = [0;1/m];