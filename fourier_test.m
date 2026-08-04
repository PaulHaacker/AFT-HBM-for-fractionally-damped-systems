function [t,x,xdot,x_proj,x_proj_dot] = fourier_test(fun,fun_derivative,omega,n,N,L)
% A test for the discrete Fourier transform
%
% Input
%  fun            : function giving x(t) as a function of t
%  fun_derivative : function giving the derivative dx/dt(t)
%  omega          : radial frequency
%  n              : number of states in x(t)
%  N              : harmonic order
%  L              : number of sample points
%
% Output
%  t              : vector of sample time instants
%  x              : vector of sample points of x
%  xdot           : vector of sample points of dx/dt
%  x_proj         : x_proj = F_N^(-1)(X), where X = F_N(x)
%  x_proj_dot     : x_proj_dot = F_N^(-1)(DX)

[V,W,D] = fourier_matrices(omega,n,N,L);
T = 2*pi/omega;
dt = T/L;
t = (0:dt:T-dt)';
x = fun(t); 
xdot = fun_derivative(t);
X = W'*x;
x_proj = V*X;
x_proj_dot = V*D*X;
