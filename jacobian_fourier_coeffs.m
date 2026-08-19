function J_realFourierStacked = jacobian_fourier_coeffs(X_sol, omega_sol, sys_jac, alpha, n, N, L)
% Computes the stacked real Fourier coefficient matrices of the time-varying
% Jacobian J(t) = df/dx evaluated along the periodic orbit X_sol.
%
% Output: J_realFourierStacked is n*(2N+1) x n, where the k-th block of n
% rows is the k-th Fourier coefficient matrix of J(t), ordered as
% (J^(0).', J^(C1).', ..., J^(CN).', J^(S1).', ..., J^(SN).')

[V, W, ~, D_alpha] = fourier_matrices(omega_sol, n, N, L, alpha);
T_period   = 2*pi / omega_sol;
t_period   = (0 : T_period/L : T_period*(1-1/L))';

x_alpha_vec = V * D_alpha * X_sol;
x_alpha     = reshape(x_alpha_vec, n, L)';

x_vec = V * X_sol;
x     = reshape(x_vec, n, L)';

dfvecdx_stacked = zeros(n*L, n);
for ii = 1:L
    [~, dfdx_i,~] = sys_jac(t_period(ii), x(ii,:)', x_alpha(ii,:)', omega_sol);
    dfvecdx_stacked(n*(ii-1)+1 : ii*n, :) = dfdx_i;
end

J_realFourierStacked = W.' * dfvecdx_stacked;
