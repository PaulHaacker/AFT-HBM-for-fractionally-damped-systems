function Lambdas = getFloquetExponents(om, X, sys_jac, alpha, n, N, L, ...
                                       real_interval, imag_interval, num_points, debug, SearchMethod, N_Hill)
% Computes Floquet exponents along a HBM continuation curve [om, X] of a fractionally damped system of the form \dot x= f(t, x, D^\alpha x) with linearization dynamics \dot y = J(t)y + C0 D^{\alpha}y, i.e. we assume f is linear in D^\alpha x.
%
% Strategy:
%   om(1)   : full search (grid search or homotopy search, see SearchMethod)
%   om(kk>=2): warm-start Newton via FractionalHillZeros_Warmstart using the
%             exponents found at om(kk-1) as initial guesses; falls back to
%             a full search (see SearchMethod) if the count changed, the
%             warm-start returned nothing, or two exponents collapsed onto
%             the same alias group.
%
% Inputs:
%   om, X          : output of arclength_continuation_HBM
%   sys_jac        : system description in the form sys_jac(t,x,x_frac,omega) with right-hand side
%                    and Jacobians (wrt x and x_frac) in dependence of omega
%   alpha          : fractional order
%   n, N, L        : state dimension, number of harmonics, sample points
%   real_interval  : [re_min, re_max] for initial grid search (default [-20, 0])
%   imag_interval  : [im_min, im_max] for initial grid search (default [-1/2, 1/2]*om(1))
%   num_points     : grid resolution per dimension for initial search (default 5)
%   debug          : pass 'debugging' to enable diagnostic prints and plots
%   SearchMethod   : 'HomotopySearch' (default) or 'GridSearch' — method used
%                    for the full search at om(1) and for every fallback
%                    search during the continuation
%   N_hill         : Fourier truncation for Hill's method (default empty, resulting in 
%                    taking just the number of harmonics N from the HBM solution)
%
% Output:
%   Lambdas : cell array, Lambdas{kk} is a column vector of Floquet exponents at om(kk)

if nargin < 13 ,  N_Hill  = [];      end
if N_Hill < N,  
    N_Hill  = [];
    warning('N_Hill must be at least as large as N, continue with N = N_Hill')
end
if nargin < 12 || isempty(SearchMethod),  SearchMethod  = 'HomotopySearch';      end
if nargin < 11 || isempty(debug),         debug         = '';                    end
if nargin < 10 || isempty(num_points),    num_points    = 5;                     end
if nargin < 9  || isempty(imag_interval), imag_interval = [-1/2, 1/2]*om(1);     end
if nargin < 8  || isempty(real_interval), real_interval = [-20, 0];              end
do_debug = strcmpi(debug, 'debugging');

addpath(fullfile(fileparts(mfilename('fullpath')), 'Hill for fractionally damped systems'))
num_om  = length(om);
Lambdas = cell(num_om, 1);

% extract C0 for homotopy search 
[~,~,C0] = sys_jac(0, zeros(n,1), zeros(n,1), om(1)); % using that right-hand side is linear in D^\alpha x, so C0 is constant

% --- kk = 1: full search ---
fprintf('getFloquetExponents: %s at om(1) = %.4f ... ', SearchMethod, om(1))
% tic
Lambdas{1} = run_full_search(SearchMethod, om(1), X(:,1), sys_jac, alpha, n, N, L, ...
                              C0, real_interval, imag_interval, num_points, N_Hill);
fprintf('%.2f s, %d exponents found\n', toc, length(Lambdas{1}))

% --- kk >= 2: warm-start Newton ---
tol_group_dup = 1e-4; % tolerance (in lambda-space) below which two exponents are considered aliases of the same group
n_exp_ref = length(Lambdas{1});
for kk = 2:num_om
    Mat_norm    = build_mat_norm(om(kk), X(:,kk), sys_jac, alpha, n, N, L, N_Hill);
    Lambdas{kk} = HillZeros_Warmstart(Mat_norm, Lambdas{kk-1});
    
    % rerun homotopy search if count changed, warm-start returned nothing,
    % or two of the warm-started exponents collapsed onto the same group
    % (the empty case must be caught explicitly: 0==0 would otherwise skip the fallback)

    has_dup = has_alias_duplicate(Lambdas{kk}, om(kk), tol_group_dup);
    if isempty(Lambdas{kk}) || length(Lambdas{kk}) ~= n_exp_ref || has_dup
        if has_dup && length(Lambdas{kk}) == n_exp_ref
            fprintf('  kk=%d: two exponents collapsed onto the same group, rerunning %s ... ', kk, SearchMethod)
        else
            fprintf('  kk=%d: exponent count changed (%d -> %d), rerunning %s ... ', ...
                    kk, n_exp_ref, length(Lambdas{kk}), SearchMethod)
        end
        tic
        Lambdas{kk} = run_full_search(SearchMethod, om(kk), X(:,kk), sys_jac, alpha, n, N, L, ...
                                       C0, real_interval, imag_interval, num_points, N_Hill);
        fprintf('%.2f s, %d exponents found\n', toc, length(Lambdas{kk}))
    end

    % % shift all exponents into the principal strip (-om/2, om/2] 
    % shifted_imag = mod(imag(Lambdas{kk}), om(kk));
    % too_high = shifted_imag > om(kk)/2;
    % shifted_imag(too_high) = shifted_imag(too_high) - om(kk);
    % Lambdas_shifted = Lambdas{kk} - imag(Lambdas{kk})*1i  + shifted_imag*1i;

    % % test whether shifted exponents are roots of the Hill matrix
    % for kk2 = 1:length(Lambdas_shifted)
    %     lambda_test = Lambdas_shifted(kk2);
    %     det_H = det(Mat_norm(lambda_test));
    %     if abs(det_H) > 1e-13
    %         warning('getFloquetExponents: shifted exponent %d at om(kk=%d) is not a root of the Hill matrix (|det(H)| = %.3e), re-running warmstart search', ...
    %                 kk2, kk, abs(det_H))
            
    %         Lambdas_shifted(kk2) = HillZeros_Warmstart(Mat_norm, Lambdas_shifted(kk2));

    %         % J_rfs   = jacobian_fourier_coeffs(X(:,kk), om(kk), sys_jac, alpha, n, N, L);
    %         % J_cell  = real2complex_fourier_jacobian(J_rfs);
    %         % lambdas_homotopy = HillZeros_Homotopy(om(kk), alpha, J_cell,C0);

    %         % Lambdas_shifted(kk2) = lambdas_homotopy;
    %     end

    % end

    % Lambdas{kk} = Lambdas_shifted;    

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
function lambdas = run_full_search(method, omega, X_col, sys_jac, alpha, n, N, L, ...
                                    C0, real_interval, imag_interval, num_points, N_Hill)
    % dispatches to the grid search or the homotopy search for a full
    % (non-warm-started) search of the Floquet exponents at a single omega
    switch lower(method)
        case 'gridsearch'
            Mat_norm = build_mat_norm(omega, X_col, sys_jac, alpha, n, N, L, N_Hill);
            lambdas_all = HillZeros(Mat_norm, omega, real_interval, imag_interval, num_points, 'newtoncomplexscalar');
            lambdas = lambdas_all(abs(imag(lambdas_all)) <= omega/2);
        case 'homotopysearch'
            J_rfs   = jacobian_fourier_coeffs(X_col, omega, sys_jac, alpha, n, N, L);
            J_cell  = real2complex_fourier_jacobian(J_rfs);

            if ~isempty(N_Hill)
                zero_block = zeros(n);
                J_added = repmat({zero_block}, 1, 2*(N_Hill-N));  % initialise with zeros
                J_cell = [J_added, J_cell, J_added]; % pad with zeros to get 4N_Hill+1 blocks
            end

            lambdas = HillZeros_Homotopy(omega, alpha, J_cell, C0);
        otherwise
            error('getFloquetExponents:unknownSearchMethod', ...
                  'Unknown SearchMethod ''%s'' (use ''HomotopySearch'' or ''GridSearch'')', method)
    end
end

% -------------------------------------------------------------------------
function tf = has_alias_duplicate(lambdas, om, tol)
    % returns true if two entries of lambdas are aliases of the same
    % underlying Floquet exponent, i.e. they differ by (approximately) an
    % integer multiple of i*om (see HillZeros.m for the same distance
    % metric used during grid-search deflation)
    tf = false;
    for ii = 1:length(lambdas)-1
        for jj = ii+1:length(lambdas)
            k    = round((imag(lambdas(ii)) - imag(lambdas(jj))) / om);
            dist = abs(lambdas(ii) - lambdas(jj) - 1i*k*om);
            if dist < tol
                tf = true;
                return
            end
        end
    end
end

% -------------------------------------------------------------------------
function Mat_norm = build_mat_norm(omega, X_col, sys_jac, alpha, n, N, L, N_Hill)
    J_rfs    = jacobian_fourier_coeffs(X_col, omega, sys_jac, alpha, n, N, L);
    J_cell   = real2complex_fourier_jacobian(J_rfs); %   J_cell : cell array of 4N+1 complex n x n matrices ordered as {J_{-2N}, ..., J_0, ..., J_{+2N}} where J_k = 0 for N < |k| <= 2N

    [~,~,C0] = sys_jac(0, zeros(n,1), zeros(n,1), omega); % using that right-hand side is linear in D^\alpha x, so C0 is constant

    if ~isempty(N_Hill)
        zero_block = zeros(n);
        J_added = repmat({zero_block}, 1, 2*(N_Hill-N));  % initialise with zeros
        J_cell = [J_added, J_cell, J_added]; % pad with zeros to get 4N_Hill+1 blocks
    end
        
    Mat_H    = giveHill(omega, alpha, J_cell, C0);
    Mat_norm = normalizeMatrix(Mat_H, n, alpha);
end

