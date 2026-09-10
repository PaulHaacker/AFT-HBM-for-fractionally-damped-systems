%% Task 3.4
clear; clc; 
% close all
%% Parameters of [PadovanSawicki1998]
% sigma = 1;
% beta = 0.01;
% gamma = 1;
% delta = 0.15;
% alpha_vec = [ 1/2, 4/5, 1, 5/4, 3/2];
% omega_start = 0.1;
% omega_end = 2;

% L = 100;
% N = 5;
% n = 2;

%% Parameters of Workshop
% sigma = 1;
% beta = 0.04;
% gamma = 1;
% delta = 0.1;
% alpha_vec = 1;

% omega_start = 0.1;
% omega_end = 2;

% L = 100;
% N = 5;
% n = 2;

%% random parameters
sigma = 1;
beta = 0.04;
gamma = 1;
delta = 0.1;
alpha_vec = [.2, .4, .6, .8, 1];

omega_start = 0.1;
omega_end = 5;

L = 100;
N = 5;
n = 2;

%% Arclength continuation with HBM - varying alpha
figure; h = gcf; h.Name = 'Task 3.4'; hold on
for k = 1:length(alpha_vec)
    alpha = alpha_vec(k);
    X0 = zeros(n*(2*N+1),1);

    [om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1);

    A = zeros(length(om),1);
    for i=1:length(om)
      [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
      A(i) = max(x(:,1));
    end
    plot(om,A,'-','DisplayName',sprintf('$\\alpha = %g$',alpha))
end
hold off
xlabel('$\omega$','Interpreter','latex')
ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
title('Nonlinear FRF of the Duffing system for different $\alpha$','Interpreter','latex')
legend('Interpreter','latex')

%% Arclength continuation with HBM - varying delta
sigma = 1;
beta = 0.04;
gamma = 1;
delta_vec = [0.1, 0.09, 0.07, 0.05,0.0001,0];
alpha = 1;

omega_start = 0.1;
omega_end = 2;

L = 100;
N = 5;
n = 2;

figure; h = gcf; h.Name = 'Task 3.4'; hold on
for k = 1:length(delta_vec)
    delta = delta_vec(k);
    X0 = zeros(n*(2*N+1),1);

    [om,X] = arclength_continuation_HBM(@(t,x,x_frac,omega) duffing_jac_frac(t,x,x_frac,sigma,beta,gamma,delta,omega),alpha,omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1);

    A = zeros(length(om),1);
    for i=1:length(om)
      [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
      A(i) = max(x(:,1));
    end
    plot(om,A,'-','DisplayName',sprintf('$\\delta = %g$',delta))
end
hold off
xlabel('$\omega$','Interpreter','latex')
ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
title(['Nonlinear FRF of the Duffing system for different $\delta$ and $\alpha = $', num2str(alpha)],'Interpreter','latex')
legend('Interpreter','latex')