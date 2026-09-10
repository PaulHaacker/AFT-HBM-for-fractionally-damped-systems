%% MMS_standardDuffingForced.m
% Multiple-scales (MMS) analysis of the standard forced Duffing oscillator.
%
% system equations:
%   \ddot x + \epsilon \beta \dot x + \omega_0^2 x + \epsilon \mu_{nl} x^3 = \epsilon f \cos(\Omega t)
%
% 0-order approximation: x(t) = a(t) cos(\Omega t - \gamma(t)), with slow dynamics of a(t), gamma(t) governed by the amplitude/phase equations below.
% Amplitude/phase equations:
%   da/dt    = -beta/2 * a + f/(2*omega0) * sin(gamma)
%   dgamma/dt = sigma - 3/8 * (mu_nl/omega0) * a^2 + f/(2*omega0*a) * cos(gamma)
%
% Equilibria (a,sigma) satisfy
%   0 = a^2 * ( (beta/2)^2 + (sigma - 3/8*(mu_nl/omega0)*a^2)^2 ) - (f/(2*omega0))^2
%
% This equation is traced with pseudo-arclength continuation (see
% arclength_continuation.m) in the detuning parameter sigma. Note that
% arclength_continuation.m rejects corrector steps with mu <= 0 (it was
% written for a continuation parameter that must stay positive, e.g. an
% excitation frequency). Since sigma is allowed to be negative here, the
% continuation is run on a shifted parameter eta = sigma + offset that
% stays positive over [sigma_start, sigma_end], and sigma is recovered
% afterwards as eta - offset.
%
% Stability of each equilibrium is obtained from the eigenvalues of the
% Jacobian of the (a,gamma) dynamics, evaluated using sin(gamma), cos(gamma)
% expressed algebraically in terms of (a,sigma) via the equilibrium
% conditions (no need to track gamma explicitly).

clear; clc; close all

%% Parameters (standard hardening Duffing oscillator, weak damping/forcing)
beta   = 0.1;   % damping coefficient
f      = 0.2;   % forcing amplitude
omega0 = 1;     % natural frequency
mu_nl  = 1;     % cubic nonlinearity coefficient (hardening for mu_nl > 0)

sigma_start = -3;
sigma_end   = 3;

%% Initial guess: low-amplitude solution at sigma_start
% F(a,sigma) is, in terms of b = a^2, a cubic k^2*b^3 - 2*sigma*k*b^2 +
% (sigma^2+(beta/2)^2)*b - c^2 = 0 with k = 3/8*mu_nl/omega0, c = f/(2*omega0).
k = 3/8 * mu_nl/omega0;
c = f/(2*omega0);
poly_coeffs = [k^2, -2*sigma_start*k, sigma_start^2 + (beta/2)^2, -c^2];
b_roots = roots(poly_coeffs);
b_roots = b_roots(abs(imag(b_roots)) < 1e-9 & real(b_roots) > 0);
a0 = sqrt(min(real(b_roots)));

%% Pseudo-arclength continuation of the equilibria curve
% continue in eta = sigma + offset, kept positive throughout, since
% arclength_continuation.m rejects corrector steps with mu <= 0
tol         = 1e-10;
maxiter     = 100;
step_start  = 0.01;
step_max    = 0.05;

offset    = max(0, 1 - sigma_start);
eta_start = sigma_start + offset;
eta_end   = sigma_end   + offset;

[eta_vec, a_vec] = arclength_continuation(...
    @(x,eta) MMS_equilibrium(x, eta - offset, beta, f, omega0, mu_nl), ...
    eta_start, eta_end, a0, tol, maxiter, step_start, step_max);
sigma_vec = eta_vec - offset;

%% Stability via eigenvalues of the (a,gamma) Jacobian
num_pts       = length(sigma_vec);
is_stable_MMS = false(1, num_pts);
for kk = 1:num_pts
    a   = a_vec(kk);
    sig = sigma_vec(kk);
    J   = MMS_jacobian(a, sig, beta, f, omega0, mu_nl);
    is_stable_MMS(kk) = all(real(eig(J)) < 0);
end

%% Do full HBM continuation for comparison

% convert parameters

epsilon = 0.01; % small parameter for weak damping/forcing/nonlinearity

sigma = omega0^2; % stiffness coefficient
delta = epsilon*beta; % damping coefficient
beta = epsilon*mu_nl; % nonlinear stiffness coefficient
gamma = epsilon*f; % forcing amplitude
alpha = 1; % fractional order of the derivative - here regular damping, so alpha = 1

omega_start = omega0 + sigma_start*epsilon; % start frequency for HBM continuation
omega_end = omega0 + sigma_end*epsilon; % end frequency for HBM continuation

N = 10;

N_Hill = [];

L = 20*floor((2*N+1)*1.2); % number of sample points of one period
n = 2;

%% Arclength continuation with HBM
X0 = zeros(n*(2*N+1),1);

[om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-3,1e-2);

A_HBM = zeros(length(om),1);
for i=1:length(om)
  [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
  A_HBM(i) = max(x(:,1));
end

% figure
% plot(om,A_HBM,'-','DisplayName',sprintf('$\\alpha = %g$',alpha))

% hold off
% xlabel('excitation frequency $\omega$','Interpreter','latex')
% ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
% title('Nonlinear FRF of the Duffing system','Interpreter','latex')
% legend('Interpreter','latex')

%% computing Floquet exponents along continuation
real_interval = [-20, 0];
imag_interval = [-1/2, 1/2]*om(1);
num_points = 5; % grid resolution for initial search of Floquet exponents at om(1) via FractionalHillZeros
Lambdas = getFloquetExponents(om, X, @(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),...
                                alpha, n, N, L, ...
                                real_interval, imag_interval, num_points, [], [], N_Hill);

%% simultaneous plotting
% classify stability from Floquet exponents
is_stable_hbm   = false(length(om), 1);
is_unstable_hbm = false(length(om), 1);
for kk = 1:length(om)
    if ~isempty(Lambdas{kk})
        if all(real(Lambdas{kk}) < 0)
            is_stable_hbm(kk)   = true;
        else
            is_unstable_hbm(kk) = true;
        end
    end
end
is_unknown_hbm = ~is_stable_hbm & ~is_unstable_hbm;

% compute excitation frequency for mms from deturning parameter sigma
Om_mms = omega0 + epsilon*sigma_vec;

% % plot: interactive slider view of the continuation curve and the
% % Floquet exponents in the complex plane at the selected point
% plotFloquetInteractive(om, A_HBM, Lambdas, is_stable_hbm, is_unstable_hbm, is_unknown_hbm,['Nonlinear FRF of duffing system for $\alpha = $', num2str(alpha)])

%% Combined plot: HBM vs. MMS, one window
% line style encodes stability (solid = stable, dashed = unstable),
% color encodes the method (HBM vs. MMS)
color_hbm = [0, 0.4470, 0.7410];    % MATLAB default blue
color_MMS = [0.8500, 0.3250, 0.0980]; % MATLAB default orange

figure; hold on
plot_stability_curve(om, A_HBM, is_stable_hbm, color_hbm)
plot_stability_curve(Om_mms, a_vec, is_stable_MMS, color_MMS)

% proxy handles so the legend shows one clean entry per method
h_hbm = plot(nan, nan, '-', 'Color', color_hbm, 'LineWidth', 1.5, 'DisplayName', 'HBM');
h_MMS = plot(nan, nan, '-', 'Color', color_MMS, 'LineWidth', 1.5, 'DisplayName', 'MMS');
legend([h_hbm, h_MMS], 'Location', 'best', 'Interpreter', 'latex')
hold off

xlabel('$\Omega$', 'Interpreter', 'latex')
ylabel('$\max\, x(t)$', 'Interpreter', 'latex')
title(['Forced Duffing oscillator: $\ddot x + \varepsilon \beta \dot x + \omega_0^2 x + \varepsilon \mu_{nl} x^3 = \varepsilon f \cos(\Omega t)$, $\varepsilon = $ ', num2str(epsilon)], 'Interpreter', 'latex')

%% -------------------------------------------------------------------
function [F, dFdx, dFdmu] = MMS_equilibrium(x, mu, beta, f, omega0, mu_nl)
% Residual and Jacobians of the MMS equilibrium equation for use with
% arclength_continuation.m. x = a (amplitude), mu = sigma (detuning).
a   = x(1);
sig = mu;
k   = 3/8 * mu_nl/omega0;
c   = f/(2*omega0);
D   = sig - k*a^2;

F     = a^2*((beta/2)^2 + D^2) - c^2;
dFdx  = 2*a*((beta/2)^2 + D^2) - 4*k*a^3*D;
dFdmu = 2*a^2*D;
end

%% -------------------------------------------------------------------
function J = MMS_jacobian(a, sig, beta, f, omega0, mu_nl)
% Jacobian of the (a,gamma) dynamics at an equilibrium (a,sig), with
% sin(gamma), cos(gamma) eliminated algebraically via the equilibrium
% conditions:
%   sin(gamma) =  beta*a*omega0/f
%   cos(gamma) = -2*omega0*a*D/f,   D = sig - 3/8*(mu_nl/omega0)*a^2
k = 3/8 * mu_nl/omega0;
D = sig - k*a^2;

J = [ -beta/2,        -a*D;
       D/a - 2*k*a,   -beta/2 ];
end

%% -------------------------------------------------------------------
function plot_stability_curve(x, y, is_stable, color)
% Plots (x,y) as contiguous solid segments where is_stable is true and
% dashed segments where it is false, all in the given color. Segments
% are excluded from the legend (HandleVisibility off); add a proxy plot
% for the legend entry instead.
x = x(:); y = y(:); is_stable = logical(is_stable(:));
n = length(x);
seg_start = 1;
for kk = 2:n+1
    if kk > n || is_stable(kk) ~= is_stable(seg_start)
        idx = seg_start:min(kk, n); % include one extra point for a continuous line
        if is_stable(seg_start)
            line_style = '-';
        else
            line_style = '--';
        end
        plot(x(idx), y(idx), 'LineStyle', line_style, 'Color', color, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off')
        seg_start = kk;
    end
end
end
