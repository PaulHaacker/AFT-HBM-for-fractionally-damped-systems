function  [V,W,D,D_alpha] = fourier_matrices(omega,n,N,L,alpha);
% matrices for real DFT and IDFT
% Remco Leine, INM, 2023
% inputs:
% omega ... frequency
% n ... number of states
% N ... number of harmonics
% L ... number of sample points
% alpha ... order of fractional damping (optional, only needed for D_alpha)
% outputs:
% V ... matrix for real DFT
% W ... matrix for real IDFT
% D ... matrix for derivative in Fourier space
% D_alpha ... matrix for fractional derivative in Fourier space

if nargin<5
    alpha = 1;
end

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

Omega_alpha = kron(diag(omegavec.^alpha),I);
c_alpha = cos(alpha*pi/2);
s_alpha = sin(alpha*pi/2);
D_alpha = [zeros(n,n) zeros(n,n*N) zeros(n,n*N);
     zeros(n*N,n) c_alpha*Omega_alpha s_alpha*Omega_alpha;
     zeros(n*N,n) -s_alpha*Omega_alpha c_alpha*Omega_alpha];


