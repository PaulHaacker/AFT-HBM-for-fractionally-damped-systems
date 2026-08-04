function [t,x,Phi_t1] = simulate_x_Phi(sys_jac,tspan,x0,odeopts)
% Simulation of the solution with computation of the fundamental solution 
% matrix by integrating the variational equation
%
% Input
%  sys_jac  : a system description giving right-hand side and Jacobian
%  tspan    : [t0 t1]
%  x0       : initial point
%  odeopts  : options for ODE solver
%
% Output
%  t        : time vector
%  x        : state vector in R^n, x(i,1:n) is the state at time t(i)
%  Phi_t1   : fundamental solution matrix at time t = t1

n = length(x0);
Phi0 = eye(n);

z0 = [x0;mat2vec(Phi0)];

[t,z] = ode45(@(t,z) sys_ext(t,z,sys_jac),tspan,z0,odeopts);

x = z(:,1:n);
Phi_t1 = vec2mat(z(end,n+1:end));
    

