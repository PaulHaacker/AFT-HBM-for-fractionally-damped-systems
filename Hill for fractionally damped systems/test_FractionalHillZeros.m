%% script testing the function FractionalHillZeros
% considering the system D^\alpha x = A(t)x(t), where here 
% A(t) = a + b\sin(\omega t) and everything scalar

close all

%% parameter stash
% a = .4;
% b = [1, 1.5, 2, 2.5, 3];

% a= -.66;
% b = [1.65,1.74,1.75];

% a = .12;
% b = 2.51;

% a = -1; % system parameter
% b = 2.5; % system parameter
% alpha = .5;
% omega = 1;


%% test parameters
% in order to test the algorithm, use a alpha=1-st order system that is
% stable and has floquet exponents lambda = a + ki, where k \in Z

a = -1; % system parameter
b = 2.5; % system parameter
alpha = 1; % DO NOT CHANGE - for testing purposes, a linear EVP is considered
omega = 1;

%% calling FractionalHillZeros

N = 10; % (2N+1) is the number of fourier coefficients considered
x_interval = [-10, 10]; % real interval for grid search
y_interval = [-10, 10]; % imaginary interval for grid search
num_points = 50; % Number of grid points to use for the search in each dimension.

% Fourier coefficients are available analytically
J = cell(1,4*N+1); % create empty cell
for kk = 1:length(J)
    J{kk} = 0;
end
J{2*N} = 1i*b/2; % J_{-1}
J{2*N+1} = a; % J_0
J{2*N+2} = -1i*b/2; % J_1

Mat_H_N_alpha = giveFracHill(omega, alpha, J);
Mat_norm = normalizeMatrix(Mat_H_N_alpha,length(J{1}),alpha);
fprintf('running FractionalHillZeros with method newtoncomplexscalar \n')
tic
lambdas_complexscalar = FractionalHillZeros(Mat_norm, x_interval, y_interval, num_points,'newtoncomplexscalar');
toc
fprintf('found a nummber of %.0f eigenvalues \n',length(lambdas_complexscalar))
fprintf('running FractionalHillZeros with method newtoncomplexscalarvectorized \n')
tic
lambdas = FractionalHillZeros(Mat_norm, x_interval, y_interval, num_points,'newtoncomplexscalarvectorized');
toc
fprintf('found a nummber of %.0f eigenvalues \n',length(lambdas))
% the returned value lambda is the numerically found lambda.

%% test plot
% Real and imaginary grids
real_grid = linspace(x_interval(1), x_interval(2), 1000);
imag_grid = linspace(y_interval(1), y_interval(2), 1000);
[real_mesh, imag_mesh] = meshgrid(real_grid, imag_grid);

% Compute |det M(lambda)|
absDet_mesh = zeros(size(real_mesh));
for kk = 1:length(real_grid)
    for jj = 1:length(imag_grid)
        absDet_mesh(jj, kk) = abs(det(Mat_H_N_alpha(real_grid(kk) + 1i * imag_grid(jj))));
    end
end

% Plot
figure
surf(real_mesh, imag_mesh, log(absDet_mesh+1), 'EdgeColor', 'none','FaceAlpha',0.8)
colorbar
set(gca,'ColorScale','log')
xlabel('real($\lambda$)')
ylabel('imag($\lambda$)')
title({'determinant of the fractional hill matrix $\log(|\det H_N^{(\alpha)}(\lambda)|+1)$,',...
    ' System ${}^CD_{-\infty}^\alpha x = (a+b\sin(\omega  t))x$',...
    ['Parameter: $a = ', num2str(a), '$, $b = ', num2str(b),...
    '$, $\alpha = ', num2str(alpha), '$, $\omega = ', num2str(omega),...
    '$, $N = ', num2str(N), '$']})
hold on
help1 = zeros(length(lambdas),1);
for kk = length(lambdas)
    help1(kk) = log(abs(det(Mat_H_N_alpha(lambdas(kk))))+1);
end
scatter3(real(lambdas), imag(lambdas), help1, 'ro', 'filled')
eigvals = eig(Mat_H_N_alpha(0));
scatter3(real(eigvals),imag(eigvals),(zeros(size(eigvals))),'bx')