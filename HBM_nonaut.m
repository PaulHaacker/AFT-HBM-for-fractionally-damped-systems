function [X,r,conv,k] = HBM_nonaut(sys_jac,X0,omega,n,N,L,tol,maxiter)
% Harmonic Balance Method for nonautonomous systems
%
% Input
%  sys_jac     : system description in the form sys_jac(t,x)
%  X0          : initial guess for the Fourier coefficients X
%  omega       : radial frequency
%  n           : number of states
%  N           : hamonic order
%  L           : number of sample points
%  tol         : tolerance for the Newton method
%  maxiter     : maximum number of iterations for the Newton method
%
% Output
%  X           : Fourier coefficients of the periodic solution
%  r           : residue of the harmonic balance method
%  conv        : boolean vector indicating convergence
%  k           : number of iterations of the Newton method

% Remco Leine, INM, 2023
T = 2*pi/omega;
dt = T/L;
t = (0:dt:T-dt)';

[V,W,D] = fourier_matrices(omega,n,N,L);
[X,r,drdx,conv,k] = newton(@(X) residue(X),X0,tol,maxiter);

    function [r,drdx] = residue(X)
    x_vec = V*X;
    x = reshape(x_vec,n,L)';
    f_vec = zeros(n*L,1);
    dfvecdxvec = zeros(n*L);
    for i=1:L
        [f_i,dfdx_i] = sys_jac(t(i),x(i,:));
        f_vec(n*(i-1)+1:i*n) = f_i;
        dfvecdxvec(n*(i-1)+1:i*n,n*(i-1)+1:i*n) = dfdx_i;
    end
    F = W'*f_vec;
    r = D*X - F;
    drdx = D-W'*dfvecdxvec*V;
    end

end

