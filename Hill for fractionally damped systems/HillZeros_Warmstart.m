function lambdas = HillZeros_Warmstart(Mat_norm, lambda0, tol_newton, tol_dup)
% Finds zeros of the normalized Hill matrix using warm-start Newton iterations.
% No grid search — each entry of lambda0 seeds one Newton run.
%
% Inputs:
%   Mat_norm   : normalized Hill matrix handle (from normalizeMatrix)
%   lambda0    : complex column vector of initial guesses (previous step's exponents)
%   tol_newton : Newton convergence tolerance (default 1e-13)
%   tol_dup    : tolerance for duplicate detection  (default 1e-4)

if nargin < 4, tol_dup    = 1e-4;  end
if nargin < 3, tol_newton = 1e-13; end

zero_fcn = @(lambda) det(Mat_norm(lambda));
lambdas  = [];

for k = 1:length(lambda0)
    [lambda, ~, exitflag] = newton_numericalGradient(zero_fcn, lambda0(k), tol_newton, 100);
    if exitflag && (isempty(lambdas) || all(abs(lambda - lambdas) > tol_dup))
        lambdas = [lambdas; lambda];
    end
end
end

% function [root, iter, exitflag] = newton_complex_scalar(f, x0, tol, max_iter)
% h = 1e-8;
% x = x0;
% for iter = 1:max_iter
%     fx  = f(x);
%     dfx = (f(x + h) - f(x - h)) / (2*h);
%     if abs(dfx) < eps
%         exitflag = 0; root = []; return
%     end
%     x_new = x - fx / dfx;
%     if abs(x_new - x) < tol
%         root = x_new; exitflag = 1; return
%     end
%     x = x_new;
% end
% exitflag = 0; root = [];
% end
