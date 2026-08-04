function [V,W,D] = fourier_matrices(omega,n,N,L)
% matrices for real DFT and IDFT
% Remco Leine, INM, 2023
T = 2*pi/omega;
dt = T/L;
t = (0:dt:T-dt)';
I = eye(n);
omegavec = (1:N)'*omega;
Gamma = [ones(L,1) cos(t*omegavec') sin(t*omegavec')];
V = kron(Gamma,eye(n));

W = zeros(size(V));
W(:,1:n) = V(:,1:n)/L;
W(:,n+1:end) = 2*V(:,n+1:end)/L; 

Omega = kron(diag(omegavec),I);
D = [zeros(n,n) zeros(n,n*N) zeros(n,n*N);
     zeros(n*N,n) zeros(n*N) Omega;
     zeros(n*N,n) -Omega zeros(n*N)];

