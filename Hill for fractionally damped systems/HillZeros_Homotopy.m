function lambdas = HillZeros_Homotopy(omega, alpha, J_cell, C0, max_iter, tol_newton, tol_dup)
% HillZeros_Homotopy finds Floquet exponents of a Hill
% problem of fractionally damped systems via homotopy continuation in alpha from alpha=1 to the target alpha.
%
% At alpha=1 the Floquet exponents are the eigenvalues of Mat_H_1(0)
% (principal strip -omega/2 < Im(lambda) <= omega/2).  alpha is then
% stepped to the target value; at each step Newton is warm-started from
% the previous exponents.
%
% Inputs:
%   omega      : excitation frequency
%   alpha      : target fractional order
%   J_cell     : Fourier coefficients of J(t), same format as giveHill
%   max_iter   : maximum number adaptive stepsize reduction (default 20)
%   tol_newton : Newton convergence tolerance           (default 1e-13)
%   tol_dup    : duplicate-detection tolerance          (default 1e-4)
%
% Output:
%   lambdas : column vector of principal Floquet exponents at target alpha

if nargin < 7, tol_dup    = 1e-4;  end
if nargin < 6, tol_newton = 1e-13; end
if nargin < 5, max_iter    = 40;    end

% --- alpha=1 starting point: eigenvalues of H_1(0) ---
% For alpha=1: det(H(lambda))=0  <=>  H_1(0)*v = lambda*v

% attention! we assume that the underlying system is a mechanical system and linear in the fractional derivative.

% extract system dimension
help = size(cell2mat(J_cell));
n = help(1);           % n - system size
f = n/2;           % f - number of degrees of freedom of mechanical system
N = (help(2)/n-1)/4;   % N - number of frequencies considered in truncation
if N~=int8(N)
    error('wrong dimension of J')
end

J_cell_initial = J_cell;

% use that in the case alpha =1 we just have regular damping, but we must bring the damping matrix C0 into the Fourier coefficients of J(t) 
C0_tilde = C0(f+1:end,1:f); % assuming that C0 = [0, 0; C0_tilde, 0]
J_cell_initial{2*N+1} = J_cell_initial{2*N+1} + [zeros(f,2*f); zeros(f,f), C0_tilde];

Mat_H1   = giveHill(omega, 1, J_cell_initial, zeros(size(C0)));
eig_all  = eig(Mat_H1(0));
in_strip = imag(eig_all) > -omega/2 & imag(eig_all) <= omega/2;
lambdas  = eig_all(in_strip);
fprintf('FractionalHillZeros_Homotopy: %d exponents in principal strip at alpha=1\n', length(lambdas))


% % DEBUGGING!
% if omega > 2.3
%     1;
% end

if isempty(lambdas)
    tolerance_principalStrip = 1e-2; % relaxation tolerance for principal strip check, to avoid losing exponents due to numerical error
    in_strip = imag(eig_all) > -omega/2 + tolerance_principalStrip & imag(eig_all) <= omega/2 + tolerance_principalStrip;
    lambdas = eig_all(in_strip);
    if length(lambdas)> n
        lambdas = lambdas(imag(lambdas) >= 0);
        if length(lambdas) ~= n
            warning('unexpected number of Floquet exponents in principal strip at alpha=1: %d (expected %d), returning empty vector', length(lambdas), n)
            lambdas = [];
            return
        end
    end
    fprintf('Found %d exponents in relaxed principal strip at alpha=1\n', length(lambdas))
end

if alpha == 1
    return
end

delta_beta = 0.02;
delta_beta_max = 0.1;
beta = 1;
goal = false;

while ~goal

    conv = false;
    m = 0;
    while ~conv && m < max_iter
        m = m + 1;
        beta_new = beta - delta_beta;
        if beta_new < alpha
            beta_new = alpha;
            goal = true;
        end

        Mat_H_beta_new = giveHill(omega, beta_new, J_cell, C0);
        Mat_H_beta_new = normalizeMatrix(Mat_H_beta_new,length(J_cell{1}),beta_new);
        zero_fcn = @(lam) det(Mat_H_beta_new(lam));

        lambdas_new = [];
        for k = 1:length(lambdas)
            [lam, ~, exitflag] = newton_numericalGradient(zero_fcn, lambdas(k), tol_newton, 100);
            if exitflag && (isempty(lambdas_new) || all(abs(lam - lambdas_new) > tol_dup))
                lambdas_new = [lambdas_new; lam]; %#ok<AGROW>
            end
        end

        if isempty(lambdas_new) % Newton did not converge
            % warning('FractionalHillZeros_Homotopy: all exponents lost at alpha=%.4f (step %d/%d)', ...
                    % beta_new, m, max_iter)
            delta_beta = delta_beta / 2;
            goal = false;
        else 
            lambdas = lambdas_new;
            beta = beta_new;
            conv = true;
            delta_beta = min(delta_beta * 1.5, delta_beta_max);
        end
    end

    if ~conv
        warning('FractionalHillZeros_Homotopy: no convergence after %d step reductions, returning empty vector of lambdas', ...
                max_iter)
        lambdas = [];
        return
    end
end


% % --- homotopy: alpha = 1 -> target ---
% N_alpha = 20;  
% alpha_vals = linspace(1, alpha, N_alpha + 1);

% for ka = 2:length(alpha_vals)
%     alpha_k  = alpha_vals(ka);
%     Mat_H_k  = giveFracHill(omega, alpha_k, J_cell);
%     % Mat_H_N_alpha = giveFracHill(omega, alpha, J);
%     Mat_H_k = normalizeMatrix(Mat_H_k,length(J_cell{1}),alpha);
%     zero_fcn = @(lam) det(Mat_H_k(lam));

%     lambdas_new = [];
%     for k = 1:length(lambdas)
%         [lam, ~, exitflag] = newton_numericalGradient(zero_fcn, lambdas(k), tol_newton, 100);
%         if exitflag && (isempty(lambdas_new) || all(abs(lam - lambdas_new) > tol_dup))
%             lambdas_new = [lambdas_new; lam]; %#ok<AGROW>
%         end
%     end
%     lambdas = lambdas_new;

%     if isempty(lambdas)
%         warning('FractionalHillZeros_Homotopy: all exponents lost at alpha=%.4f (step %d/%d)', ...
%                 alpha_k, ka-1, N_alpha)
%         return
%     end
% end
% fprintf('FractionalHillZeros_Homotopy: done, %d exponents at alpha=%.4f\n', length(lambdas), alpha)
