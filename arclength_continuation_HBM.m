function [om,X] = arclength_continuation_HBM(sys_jac,om_start,om_end,X0,n,N,L,tol,maxiter,sigma_start,sigma_max,plot_fun)
% Pseudo-arclength continuation with harmonic balance method for non-autonomous systems
%
% Input 
%  sys_jac     : system description in the form sys_jac(t,x,omega) with right-hand side 
%                and Jacobian in dependence of omega
%  om_start    : first value of omega
%  om_end      : last value of omega
%  X0          : initial guess for for the Fourier coefficients X
%  n           : number of states
%  N           : number of harmonics
%  L           : number of sample points
%  tol         : tolerance
%  maxiter     : maximum number of iterations of the Newton method
%  sigma_start : starting value for the stepsize sigma
%  sigma_max   : maximal value for the stepsize sigma
%  plot_fun    : a plotting function
%
% Output
%  om          : vector with omega values
%  X           : matrix in which the i-th row are the Fourier coefficients at omega(i)

% Remco Leine, INM, University of Stuttgart, 2023

[om,X] = arclength_continuation(@(X,omega) residue(X,omega),om_start,om_end,X0,tol,maxiter,sigma_start,sigma_max,plot_fun);

    function [r,drdx,drdomega] = residue(X,omega)
    T = 2*pi/omega;
    dt = T/L;
    t = (0:dt:T-dt)';
    [V,W,D] = fourier_matrices(omega,n,N,L); 
    x_vec = V*X;
    x = reshape(x_vec,n,L)';
    f_vec = zeros(n*L,1);
    dfvecdxvec = zeros(n*L);
    for i=1:L
        [f_i,dfdx_i] = sys_jac(t(i),x(i,:),omega);
        f_vec(n*(i-1)+1:i*n) = f_i;
        dfvecdxvec(n*(i-1)+1:i*n,n*(i-1)+1:i*n) = dfdx_i;
    end
    F = W'*f_vec;
    r = D*X - F;
    drdx = D-W'*dfvecdxvec*V;
    drdomega = omega\D*X;
    end

end