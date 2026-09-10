%% WizardsHat frequency response curve
clear; clc;
 close all
%% Parameters of [AmabelliEtAl2021]

% sigma = 4;
% iota = 1.75; %0.01
% beta = .35; %0.01
% gamma = 0.2;
% delta = 0.02;
% alpha = 0.9;
% omega_start = 1;
% omega_end = 3;

% N = 5;
% N_Hill = [];

% L = 20*floor((2*N+1)*1.2); % number of sample points of one period
% n = 2;

%% random Parameters 

% sigma = 4;
% iota = 1.75; %0.01
% beta = .35; %0.01
% gamma = 0.2;
% delta = 0.05;
% alpha = .75;
% omega_start = 1;
% omega_end = 3;

% N = 10;
% N_Hill = [];

% L = 20*floor((2*N+1)*1.2); % number of sample points of one period
% n = 2;


%% difficult Parameters (stability change at fold bifurcation doesnt get recognized)

sigma = 4;
iota = 1.75; %0.01
beta = .35; %0.01
gamma = 0.2;
delta = 0.1;
alpha = .25;
omega_start = 1;
omega_end = 3;

N = 10;
N_Hill = [];

L = 20*floor((2*N+1)*1.2); % number of sample points of one period
n = 2;


%% Arclength continuation with HBM
X0 = zeros(n*(2*N+1),1);

% plot_fun = @(omega, X) plot_fun_HBM(omega, X, n, N, L);

[om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) WizardsHat_jac_frac(t,x,x_frac,sigma,iota,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1); %, plot_fun);

A_HBM = zeros(length(om),1);
for i=1:length(om)
  [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
  A_HBM(i) = max(x(:,1));
end
figure
plot(om,A_HBM,'-','DisplayName',sprintf('$\\alpha = %g$',alpha))

hold off
xlabel('$\omega$','Interpreter','latex')
ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
title('Nonlinear FRF','Interpreter','latex')
legend('Interpreter','latex')

%% computing Floquet exponents along continuation
real_interval = [-20, 0];
imag_interval = [-1/2, 1/2]*om(1);
num_points = 5; % grid resolution for initial search of Floquet exponents at om(1) via FractionalHillZeros
Lambdas = getFloquetExponents(om, X, @(t,x,x_frac,omega) WizardsHat_jac_frac(t,x,x_frac,sigma,iota,beta,gamma,delta,omega),...
                                alpha, n, N, L, ...
                                real_interval, imag_interval, num_points, [], [], N_Hill);

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
plotFloquetInteractive(om, A_HBM, Lambdas, is_stable, is_unstable, is_unknown,['Nonlinear FRF of duffing system for $\alpha = $', num2str(alpha)])

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
