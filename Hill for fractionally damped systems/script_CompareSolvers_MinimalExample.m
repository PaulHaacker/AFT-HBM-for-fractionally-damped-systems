%% script applying the fractional hill method to a minimal example system, comparing different solvers
% considering the system D^\alpha x = A(t)x(t), where here 
% A(t) = a + b\sin(\omega t) and everything scalar

clear
close all

%% parameter stash

% a = -1; % system parameter
% b = 2.2; % system parameter
% alpha = .5;
% omega = 1;


% % skewed column and infinite real ones
% a = .1;%-1; % system parameter
% b = 2.5; % system parameter
% alpha = .2;
% omega = 1;
% 
a = .1;-1; % system parameter
b = 2.5;3;2.2; % system parameter
alpha = .5;
omega = 1;

% a = -1; % system parameter
% b = 2.5; % system parameter
% alpha = .5;
% omega = 1;

%% calling FractionalHillZeros

N = 10; % (2N+1) is the number of fourier coefficients considered

% % search grid in the complex plane
% x_interval = [-5, 120]; % real interval for grid search
% y_interval = [-20, 20]; % imaginary interval for grid search
% num_points = 300; % Number of grid points to use for the search in each dimension.

% x_interval = [1, 2];[-1, 1]; % real interval for grid search
% y_interval = [0, 4]; % imaginary interval for grid search
x_interval = [1, 2];[.4, .5];[0.01,.02]; % real interval for grid search
y_interval = [0,10];[.99, 1.01]; % imaginary interval for grid search
num_points = 100; % Number of grid points to use for the search in each dimension.

% x_interval = [0,.1081]; % real interval for grid search
% y_interval = [0,1]; % imaginary interval for grid search
% num_points = 2; % Number of grid points to use for the search in each dimension.

% % Fourier coefficients are available analytically
J = cell(1,4*N+1); % create empty cell
for kk = 1:length(J)
    J{kk} = 0;
end
J{2*N} = 1i*b/2; % J_{-1}
J{2*N+1} = a; % J_0
J{2*N+2} = -1i*b/2; % J_1

Mat_H_N_alpha = giveFracHill(omega, alpha, J);

%% own implementation using scalar representation det F(z)=0
tic
lambdas_manual = FractionalHillZeros(Mat_H_N_alpha, x_interval, y_interval, num_points);
time_manual = toc;

%% NLEIGS of S. Guettel, R. Van Beeumen, K. Meerbergen, and W. Michiels.
tic
NLEP.Fun = Mat_H_N_alpha;
NLEP.n = size(Mat_H_N_alpha(0),1);
Sigma = [x_interval(1);x_interval(1);x_interval(2);x_interval(2)]...
    + 1i*[y_interval(1);y_interval(2);y_interval(2);y_interval(1)];

% options.maxit = 300;
% [X,lambdas_NLEIGS,res,info] = nleigs(NLEP,Sigma,[],options);

[X,lambdas_NLEIGS,res,info] = nleigs(NLEP,Sigma,[]);

time_NLEIGS = toc;

%% plot only eigenvalues in complex plane
fig_ComplexPlane = figure;
if ~isempty(lambdas_manual)
handle_manual = scatter(real(lambdas_manual), imag(lambdas_manual), 'ro', 'filled','DisplayName',['grid search, \#$\lambda = $ ', num2str(length(lambdas_manual)), ', time = ',num2str(time_manual)]);
end
hold on
handle_NLEIGS = scatter(real(lambdas_NLEIGS), imag(lambdas_NLEIGS), 'k*','DisplayName',['NLEIGS, \#$\lambda = $ ', num2str(length(lambdas_NLEIGS)), ', time = ',num2str(time_NLEIGS)]);
% pgon = polyshape(real(Sigma),imag(Sigma));
% plot(pgon,'FaceColor','red','FaceAlpha',0.1,'DisplayName','target set')
fill(real(Sigma), imag(Sigma), 'red', 'FaceAlpha', 0.1, 'DisplayName', 'target set')
% plot(real(Sigma),imag(Sigma),'DisplayName','target set')
xlabel('real($\lambda$)')
ylabel('imag($\lambda$)')
title({'Solutions of the Fractional Hill Matrix $\det H_N^{(\alpha)}(\lambda)=0$,',...
    ' System ${}^CD_{-\infty}^\alpha x = (a+b\sin(\omega  t))x$',...
    ['Parameter: $a = ', num2str(a), '$, $b = ', num2str(b),...
    '$, $\alpha = ', num2str(alpha), '$, $\omega = ', num2str(omega),...
    '$, $N = ', num2str(N), '$'], ...
    ['Matrix size $2N+1 = $',num2str(2*N+1),', number of solutions found \#$\lambda$ ']})
legend
% save figure
% savePlot(fig_ComplexPlane,sprintf('Det_H_N_a_%.2f_b_%.2f_alpha_%.2f_omega_%.2f_N_%d_ComplexPlane', ...
%        a, b, alpha, omega, N));

%% construct the identified time-signal

% if ~isempty(lambdas)
%     [~, index] = min(abs(imag(lambdas)));
%     lambda_sol = lambdas(index);
% 
%     p_vec = null(Mat_H_N_alpha(lambda_sol));
%     omega_vec = (-N:N)';
% 
%     soln = @(t) exp(lambda_sol*t).*sum(p_vec.*exp(1i*omega_vec*t));
% 
%     t_vec = -1:0.01:10;
%     figure
%     hold on
%     plot(t_vec,real(soln(t_vec)),'DisplayName','$\mathrm{Re}(y(t))$');
%     plot(t_vec,imag(soln(t_vec)),'DisplayName','$\mathrm{Im}(y(t))$');
%     legend
%     figure
%     hold on
%     plot(t_vec,abs(soln(t_vec)),'DisplayName','$|(y(t))|$');
%     plot(t_vec,angle(soln(t_vec)),'DisplayName','$\mathrm{arg}(y(t))$');
%     legend
% end

%% functions
function savePlot(fig_handle,name)
% saves figure with handle fig_handle as matlab .fig and .png with the name
% specified as 'name' into the folder 'plots'

% example call:
% savePlot(gcf,sprintf('Det_H_N_a_%.2f_b_%.2f_alpha_%.2f_omega_%.2f_N_%d', ...
%        a, b, alpha, omega, N));

% Construct filenames
fig_filename = fullfile('plots', [name,'.fig']);
png_filename = fullfile('plots', [name,'.png']);

% Save as .fig and .png
savefig(fig_handle,fig_filename);
saveas(fig_handle, png_filename);
end