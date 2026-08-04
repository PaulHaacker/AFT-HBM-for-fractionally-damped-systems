function [t,x] = fourier_coeff_to_time_series(X,omega,n,N,L)
% fourier_coeff_to_time_series maps X to (t,x)
% such that we can plot x(t) using plot(t,x)
% 
% Input
%  X     : vector of Fourier coefficients
%  omega : radial frequency
%  n     : dimension of the signal x
%  N     : harmonic order
%  L     : number of sample points
%
% Output
%  t     : vector of sample time-instants
%  x     : matrix in which the i-th row is the state time time t_i
T = 2*pi/omega;
dt = T/L;

% create vector of time-instants
t = (0:dt:T-dt)';

% compute Fourier matrices
[V,W,D] = fourier_matrices(omega,n,N,L);

% compute vector of sample points (inverse Fourier transform)
x_vec = V*X;

% re-order the vector of sample points in a matrix, i.e. the time-series
x = reshape(x_vec,n,L)';


