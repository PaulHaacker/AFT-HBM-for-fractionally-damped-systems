function Lambdas = getFloquetExponents(om, X, sys_jac, alpha, n, N, L, ...
                                       real_interval, imag_interval, num_points, debug)
% Computes Floquet exponents along a HBM continuation curve [om, X] of a fractionally damped system of the form \dot x= f(t, x, D^\alpha x) with linearization dynamics \dot y = J(t)y + C0 D^{\alpha}y, i.e. we assume f is linear in D^\alpha x.
%
% Strategy:
%   om(1)   : full grid search via FractionalHillZeros, then keep all
%             exponents with |imag(lambda)| <= imag_interval(2).
%   om(kk>=2): warm-start Newton via FractionalHillZeros_Warmstart using the
%             exponents found at om(kk-1) as initial guesses — no grid search.
%
% Inputs:
%   om, X          : output of arclength_continuation_HBM
%   sys_jac        : system description in the form sys_jac(t,x,x_frac,omega) with right-hand side 
%                    and Jacobians (wrt x and x_frac) in dependence of omega
%   alpha          : fractional order
%   n, N, L        : state dimension, number of harmonics, sample points
%   real_interval  : [re_min, re_max] for initial grid search
%   imag_interval  : [im_min, im_max] for initial grid search
%   num_points     : grid resolution per dimension for initial search
%   debug          : pass 'debugging' to enable diagnostic prints and plots
%
% Output:
%   Lambdas : cell array, Lambdas{kk} is a column vector of Floquet exponents at om(kk)

if nargin < 11, debug = ''; end
do_debug = strcmpi(debug, 'debugging');

addpath(fullfile(fileparts(mfilename('fullpath')), 'Hill for fractionally damped systems'))
num_om  = length(om);
Lambdas = cell(num_om, 1);

% --- kk = 1: full grid search ---
fprintf('getFloquetExponents: grid search at om(1) = %.4f ... ', om(1))
Mat_norm    = build_mat_norm(om(1), X(:,1), sys_jac, alpha, n, N, L);
lambdas_all = HillZeros(Mat_norm, real_interval, imag_interval, num_points, 'newtoncomplexscalar');
Lambdas{1}  = lambdas_all(abs(imag(lambdas_all)) <= om(1)/2);
fprintf('found %d exponents\n', length(Lambdas{1}))

% % --- kk = 1: Homotopy search
% fprintf('FractionalHillZeros_Homotopy ... ')
% tic
% J_rfs   = jacobian_fourier_coeffs(X(:,1), om(1), sys_jac, alpha, n, N, L);
% J_cell  = real2complex_fourier_jacobian(J_rfs);
% lambdas_homotopy = FractionalHillZeros_Homotopy(om(1), alpha, J_cell);
% fprintf('%.2f s, %d exponents found\n', toc, length(lambdas_homotopy))
% Lambdas{1} = lambdas_homotopy;

% --- kk >= 2: warm-start Newton ---
n_exp_ref = length(Lambdas{1});
for kk = 2:num_om
    Mat_norm    = build_mat_norm(om(kk), X(:,kk), sys_jac, alpha, n, N, L);
    Lambdas{kk} = HillZeros_Warmstart(Mat_norm, Lambdas{kk-1});
    % rerun grid search if count changed or warm-start returned nothing
    % (the empty case must be caught explicitly: 0==0 would otherwise skip the fallback)
    if isempty(Lambdas{kk}) || length(Lambdas{kk}) ~= n_exp_ref
        fprintf('  kk=%d: exponent count changed (%d -> %d), rerunning grid search ... ', ...
                kk, n_exp_ref, length(Lambdas{kk}))
        lambdas_retry = HillZeros(Mat_norm, real_interval, imag_interval, num_points, 'newtoncomplexscalar');
        Lambdas{kk}   = lambdas_retry(abs(imag(lambdas_retry)) <= om(kk)/2);
        n_exp_ref     = length(Lambdas{kk});
        fprintf('found %d exponents\n', n_exp_ref)

        % tic
        % fprintf('  kk=%d: exponent count changed (%d -> %d), rerunning homotopy search ... ', ...
        %         kk, n_exp_ref, length(Lambdas{kk}))
        % J_rfs   = jacobian_fourier_coeffs(X(:,kk), om(kk), sys_jac, alpha, n, N, L);
        % J_cell  = real2complex_fourier_jacobian(J_rfs);
        % lambdas_homotopy = FractionalHillZeros_Homotopy(om(kk), alpha, J_cell);
        % Lambdas{kk} = lambdas_homotopy;
        % fprintf('%.2f s, %d exponents found\n', toc, length(lambdas_homotopy))
    end
    if mod(kk, 50) == 0
        fprintf('  processed %d / %d\n', kk, num_om)
    end
end
fprintf('getFloquetExponents: done\n')

if do_debug
    n_lambdas = cellfun(@length, Lambdas);
    if all(n_lambdas == n_lambdas(1))
        fprintf('DEBUG: all %d omega points have %d Floquet exponents\n', num_om, n_lambdas(1))
    else
        fprintf('DEBUG: number of Floquet exponents is NOT constant (min %d, max %d)\n', min(n_lambdas), max(n_lambdas))
        figure;
        plot(n_lambdas, 'k.')
        xlabel('index $k$', 'Interpreter', 'latex')
        ylabel('number of Floquet exponents', 'Interpreter', 'latex')
        title('DEBUG: Floquet exponent count along continuation', 'Interpreter', 'latex')
    end

    % plot |lambda(kk) - lambda(kk-1)| for each exponent group
    n_exp = n_lambdas(1);
    diffs = NaN(num_om-1, n_exp);
    for kk = 2:num_om
        lk  = Lambdas{kk};
        lkm = Lambdas{kk-1};
        m   = min(length(lk), length(lkm));
        diffs(kk-1, 1:m) = abs(lk(1:m) - lkm(1:m));
    end
    figure;
    semilogy(diffs)
    xlabel('index $k$', 'Interpreter', 'latex')
    ylabel('$|\lambda_k - \lambda_{k-1}|$', 'Interpreter', 'latex')
    title('DEBUG: step-to-step change of Floquet exponents', 'Interpreter', 'latex')
    legend(arrayfun(@(k) sprintf('group %d', k), 1:n_exp, 'UniformOutput', false))

    % plot trajectories of Floquet exponents in the complex plane
    % colour encodes position along continuation (blue=start, red=end)
    lk_re = []; lk_im = []; lk_c = [];
    for kk = 1:num_om
        lk = Lambdas{kk};
        if ~isempty(lk)
            lk_re = [lk_re; real(lk)]; %#ok<AGROW>
            lk_im = [lk_im; imag(lk)]; %#ok<AGROW>
            lk_c  = [lk_c;  repmat(kk, length(lk), 1)]; %#ok<AGROW>
        end
    end
    figure; hold on
    scatter(lk_re, lk_im, 10, lk_c, 'filled')
    xline(0, 'k--', 'LineWidth', 0.8)
    yline(0,  'k-',  'LineWidth', 0.5)
    colormap(parula); cb = colorbar;
    clim([1, length(om)])
    cb.Label.String = 'continuation index $k$';
    cb.Label.Interpreter = 'latex';
    xlabel('$\mathrm{Re}(\lambda)$', 'Interpreter', 'latex')
    ylabel('$\mathrm{Im}(\lambda)$', 'Interpreter', 'latex')
    title('DEBUG: Floquet exponents in complex plane', 'Interpreter', 'latex')
    axis tight; box on
end
end

% -------------------------------------------------------------------------
function Mat_norm = build_mat_norm(omega, X_col, sys_jac, alpha, n, N, L)
    J_rfs    = jacobian_fourier_coeffs(X_col, omega, sys_jac, alpha, n, N, L);
    J_cell   = real2complex_fourier_jacobian(J_rfs);
    [~,~,C0] = sys_jac(0, zeros(n,1), zeros(n,1), omega); % using that right-hand side is linear in D^\alpha x, so C0 is constant
    Mat_H    = giveHill(omega, alpha, J_cell, C0);
    Mat_norm = normalizeMatrix(Mat_H, n, alpha);
end

