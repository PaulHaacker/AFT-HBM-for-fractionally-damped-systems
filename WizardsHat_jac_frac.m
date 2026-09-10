function [f,dfdx,dfdx_frac] = WizardsHat_jac_frac(t,x,x_frac,sigma,iota,beta,gamma,delta,omega)
% Wizards Hat system with Jacobians
% all parameters are normalized
% comment: the function does not know the order of the fractional derivative. In the chosen structure of the AFT-HBM, we dont need it here.

% t ... time
% x ... state vector
% x_frac ... fractional derivative of x
% sigma ... linear stiffness coefficient
% iota ... quadratic coefficient
% beta ... cubic stiffness coefficient
% gamma ... forcing amplitude
% delta ... fractional damping coefficient

f = [x(2) ; 
     -delta*x_frac(1)-sigma*x(1)-iota*x(1)^2-beta*x(1)^3+gamma*cos(omega*t) ];

dfdx = [ 0, 1;
         -sigma-2*iota*x(1)-3*beta*x(1)^2, 0];

dfdx_frac = [0, 0;
          -delta, 0];  

