%% Task 3.4
clear; clc; close all
%% Parameters of [PadovanSawicki1998]

% sigma = 1;
% beta = 0.01; %0.01
% gamma = 1;
% delta = 0.15;
% alpha = 1/4;
% omega_start = 0.1;
% omega_end = 2;

% N = 10;

% L = 20*floor((2*N+1)*1.2); % number of sample points of one period
% n = 2;

%% Parameters of Workshop
% sigma = 1;
% beta = 0.04;
% gamma = 1;
% delta = 0.1;
% alpha = 1;

% omega_start = 0.1;
% omega_end = 2;

% L = 100;
% N = 5;
% n = 2;

%% random parameters

sigma = 1;
beta = -0.01; %0.01
gamma = .75;
delta = 0.15;
alpha = 1/2;
omega_start = 0.1;
omega_end = 2;

N = 15;

L = 20*floor((2*N+1)*1.2); % number of sample points of one period
n = 2;

%% Arclength continuation with HBM
X0 = zeros(n*(2*N+1),1);

[om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1);

A_HBM = zeros(length(om),1);
for i=1:length(om)
  [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
  A_HBM(i) = max(x(:,1));
end
plot(om,A_HBM,'-','DisplayName',sprintf('$\\alpha = %g$',alpha))

hold off
xlabel('$\omega$','Interpreter','latex')
ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
title('Nonlinear FRF of the Duffing system','Interpreter','latex')
legend('Interpreter','latex')

%% computing Floquet exponents along continuation
real_interval = [-20, 0];
imag_interval = [-1/2, 1/2]*om(1);
num_points = 5; % grid resolution for initial search of Floquet exponents at om(1) via FractionalHillZeros
Lambdas = getFloquetExponents(om, X, @(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),...
                                alpha, n, N, L, ...
                                real_interval, imag_interval, num_points);

%% plotting       
% classify stability from Floquet exponents
is_stable   = false(length(om), 1);
is_unstable = false(length(om), 1);
for kk = 1:length(om)
    if ~isempty(Lambdas{kk})
        if all(real(Lambdas{kk}) < 0)
            is_stable(kk)   = true;
        else
            is_unstable(kk) = true;
        end
    end
end
is_unknown = ~is_stable & ~is_unstable;

% plot: interactive slider view of the continuation curve and the
% Floquet exponents in the complex plane at the selected point
plotFloquetInteractive(om, A_HBM, Lambdas, is_stable, is_unstable, is_unknown)

%% Arclength continuation with HBM - varying alpha
% figure; h = gcf; h.Name = 'Task 3.4'; hold on
% for k = 1:length(alpha_vec)
%     alpha = alpha_vec(k);
%     X0 = zeros(n*(2*N+1),1);

%     [om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1);

%     A = zeros(length(om),1);
%     for i=1:length(om)
%       [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
%       A(i) = max(x(:,1));
%     end
%     plot(om,A,'-','DisplayName',sprintf('$\\alpha = %g$',alpha))
% end
% hold off
% xlabel('$\omega$','Interpreter','latex')
% ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
% title('Nonlinear FRF of the Duffing system for different $\alpha$','Interpreter','latex')
% legend('Interpreter','latex')

%% Arclength continuation with HBM - varying delta
% sigma = 1;
% beta = 0.04;
% gamma = 1;
% delta_vec = [0.1, 0.09, 0.07, 0.05,0.0001,0];
% alpha = 1;

% omega_start = 0.1;
% omega_end = 2;

% L = 100;
% N = 5;
% n = 2;

% figure; h = gcf; h.Name = 'Task 3.4'; hold on
% for k = 1:length(delta_vec)
%     delta = delta_vec(k);
%     X0 = zeros(n*(2*N+1),1);

%     [om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1);

%     A = zeros(length(om),1);
%     for i=1:length(om)
%       [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
%       A(i) = max(x(:,1));
%     end
%     plot(om,A,'-','DisplayName',sprintf('$\\delta = %g$',delta))
% end
% hold off
% xlabel('$\omega$','Interpreter','latex')
% ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
% title(['Nonlinear FRF of the Duffing system for different $\delta$ and $\alpha = $', num2str(alpha)],'Interpreter','latex')
% legend('Interpreter','latex')

%% local functions
function plotFloquetInteractive(om, A_HBM, Lambdas, is_stable, is_unstable, is_unknown)
%PLOTFLOQUETINTERACTIVE Slider view of a continuation curve and its Floquet exponents.
%   Left panel: FRF continuation curve (om vs. A_HBM) colored by stability,
%   with a marker showing the currently selected point.
%   Right panel: Floquet exponents (from Lambdas{kk}) in the complex plane
%   at the point currently selected by the slider.

n_om = length(om);

% fixed axis limits for the complex-plane panel so the view doesn't jump
% around as the slider moves
all_lambda = vertcat(Lambdas{:});
if isempty(all_lambda)
    re_lim = [-1, 1]; im_lim = [-1, 1];
else
    re_pad = 0.1*max(max(real(all_lambda))-min(real(all_lambda)), 1e-6);
    im_pad = 0.1*max(max(imag(all_lambda))-min(imag(all_lambda)), 1e-6);
    re_lim = [min(real(all_lambda))-re_pad, max(real(all_lambda))+re_pad];
    im_lim = [min(imag(all_lambda))-im_pad, max(imag(all_lambda))+im_pad];
end
if diff(re_lim)==0, re_lim = re_lim + [-1, 1]; end
if diff(im_lim)==0, im_lim = im_lim + [-1, 1]; end

fig = figure('Name', 'Floquet exponents along continuation', 'Position', [100 100 1100 550]);

ax1 = subplot(1,2,1); hold(ax1, 'on');
plot(ax1, om, A_HBM, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(ax1, om(is_stable),   A_HBM(is_stable),   'b.', 'MarkerSize', 8, 'DisplayName', 'all Floquet exponents $\lambda$ have $\mathrm{Re}(\lambda) <0$')
plot(ax1, om(is_unstable), A_HBM(is_unstable), 'r.', 'MarkerSize', 8, 'DisplayName', '$\exists$ Floquet exponent $\lambda$ with $\mathrm{Re}(\lambda) >0$')
plot(ax1, om(is_unknown),  A_HBM(is_unknown),  'k.', 'MarkerSize', 8, 'DisplayName', 'no Floquet exponents found')
marker = plot(ax1, om(1), A_HBM(1), 'go', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'selected point');
hold(ax1, 'off');
xlabel(ax1, '$\omega$', 'Interpreter', 'latex')
ylabel(ax1, '$\max\,_{t \in [0,T]} |x_1(t)|$', 'Interpreter', 'latex')
title(ax1, 'Nonlinear FRF of the Duffing system with fractional damping', 'Interpreter', 'latex')
legend(ax1, 'Interpreter', 'latex', 'Location', 'best')
set(ax1, 'FontSize', 12)

ax2 = subplot(1,2,2); hold(ax2, 'on');
scatter(ax2, real(all_lambda), imag(all_lambda), 10, [0.7, 0.7, 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off')
xline(ax2, 0, 'k--', 'HandleVisibility', 'off');
sc = scatter(ax2, real(Lambdas{1}), imag(Lambdas{1}), 40, 'filled');
hold(ax2, 'off');
xlim(ax2, re_lim); ylim(ax2, im_lim);
xlabel(ax2, '$\mathrm{Re}(\lambda)$', 'Interpreter', 'latex')
ylabel(ax2, '$\mathrm{Im}(\lambda)$', 'Interpreter', 'latex')
title(ax2, sprintf('Floquet exponents at $\\omega = %.4f$ ($k = 1$)', om(1)), 'Interpreter', 'latex')
set(ax2, 'FontSize', 12)
grid(ax2, 'on')

step = 1/max(n_om-1, 1);
uicontrol('Parent', fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.15, 0.02, 0.7, 0.03], 'Min', 1, 'Max', n_om, 'Value', 1, ...
    'SliderStep', [step, min(10*step, 1)], ...
    'Callback', @(src, ~) update_view(round(get(src, 'Value'))));

    function update_view(kk)
        kk  = max(1, min(n_om, kk));
        lam = Lambdas{kk};
        set(marker, 'XData', om(kk), 'YData', A_HBM(kk))
        set(sc, 'XData', real(lam), 'YData', imag(lam))
        title(ax2, sprintf('Floquet exponents at $\\omega = %.4f$ ($k = %d$)', om(kk), kk), 'Interpreter', 'latex')
    end
end