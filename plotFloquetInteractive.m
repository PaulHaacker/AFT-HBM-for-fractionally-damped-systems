function plotFloquetInteractive(om, A_HBM, Lambdas, is_stable, is_unstable, is_unknown, plot_title)
%PLOTFLOQUETINTERACTIVE Slider view of a continuation curve and its Floquet exponents.
%   Left panel: FRF continuation curve (om vs. A_HBM) colored by stability,
%   with a marker showing the currently selected point.
%   Right panel: Floquet exponents (from Lambdas{kk}) in the complex plane
%   at the point currently selected by the slider.
% plot_title: optional title for the left panel (default: 'Nonlinear FRF of the Duffing system with fractional damping')
if nargin < 7
    plot_title = 'Nonlinear FRF of the Duffing system with fractional damping';
end

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
title(ax1, plot_title, 'Interpreter', 'latex')
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