%% script applying the fractional hill method to a minimal example system
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
% alpha = .5;
% omega = 1;

a = -1;%-1; % system parameter
b = 2.5; % system parameter
alpha = .5;
omega = 1;
n = 1;

% a = .1;%-1; % system parameter
% b = 1.5; % system parameter
% alpha = .2;
% omega = 1;
% n = 1;

% a = .1;%-1; % system parameter
% b = 1.4; % system parameter
% alpha = .2;
% omega = 1;
% n = 1;

%% calling FractionalHillZeros

N = 20; % (2N+1) is the number of fourier coefficients considered

% % search grid in the complex plane
% x_interval = [-5, 120]; % real interval for grid search
% y_interval = [-20, 20]; % imaginary interval for grid search
% num_points = 300; % Number of grid points to use for the search in each dimension.

% x_interval = [0, .2];[-1, 1]; % real interval for grid search
x_interval = [-1, 1]; % real interval for grid search
y_interval = [-5, 5]; % imaginary interval for grid search
num_points = 20; % Number of grid points to use for the search in each dimension.

% x_interval = [0,.1081]; % real interval for grid search
% y_interval = [0,1]; % imaginary interval for grid search
% num_points = 2; % Number of grid points to use for the search in each dimension.

% % Fourier coefficients are available analytically
% J = cell(1,4*N+1); % create empty cell
J = repmat({zeros(n,n)},1,4*N+1);
J{2*N} = 1i*b/2; % J_{-1}
J{2*N+1} = a; % J_0
J{2*N+2} = -1i*b/2; % J_1

fprintf('Finding Eigenvalues of Fractional Hill of N=%.0f frequencies. \n',N)
tic
Mat_H_N_alpha = giveFracHill(omega, alpha, J);
Mat_norm = normalizeMatrix(Mat_H_N_alpha,length(J{1}),alpha);
lambdas = FractionalHillZeros(Mat_norm, x_interval, y_interval, num_points);
toc
fprintf('found a number of %.0f eigenvalues\n',length(lambdas))

% DEBUGGING _______________________________________________________________
% Mat_H_N_alpha = giveFracHill(omega, alpha, J);
% NLEP.Fun = Mat_H_N_alpha;
% NLEP.n = size(Mat_H_N_alpha(0),1);
% Sigma = [x_interval(1);x_interval(1);x_interval(2);x_interval(2)]...
%     + 1i*[y_interval(1);y_interval(2);y_interval(2);y_interval(1)];
% 
% % options.maxit = 300;
% % [X,lambdas_NLEIGS,res,info] = nleigs(NLEP,Sigma,[],options);
% 
% [X,lambdas,res,info] = nleigs(NLEP,Sigma,[]);

% DEBUGGING END ___________________________________________________________

%% numerical eigenvalue counting
% % gam =5i;
% gam =0.1;
% rad = .5; nc = 1000; % centre, radius, nr of nodes on circle
% h = 1e-7; % for finite diff
% f = @(z) det(Mat_norm(z));
% df = @(z) finiteDiff(f,z,h);
% 
% nbr_evs = EstNumEV(f,df,gam,rad,nc);
% fprintf('Circle with center %.2f, radius %.2f has an estimate of %.4f eigenvalues.\n', gam, rad, nbr_evs);

%% test plot of the found lambda
% Real and imaginary grids
real_grid = linspace(x_interval(1), x_interval(2), 100);
imag_grid = linspace(y_interval(1), y_interval(2), 100);
[real_mesh, imag_mesh] = meshgrid(real_grid, imag_grid);

% Compute |det M(lambda)|
absDet_mesh = zeros(size(real_mesh));
for kk = 1:length(real_grid)
    for jj = 1:length(imag_grid)
        absDet_mesh(jj, kk) = abs(det(Mat_H_N_alpha(real_grid(kk) + 1i * imag_grid(jj))));
    end
end

%% Plot
fig_surf = figure('units','normalized','outerposition',[0 0 1 1]);
set(fig_surf, 'DefaultTextFontSize', 14);
set(fig_surf, 'DefaultAxesFontSize', 14);

surf(real_mesh, imag_mesh, log(absDet_mesh), 'EdgeColor', 'none','FaceAlpha',0.8)
colorbar
set(gca,'ColorScale','log')
xlabel('real($\lambda$)')
ylabel('imag($\lambda$)')
title({'determinant of the fractional hill matrix $\log\det H_N^{(\alpha)}(\lambda)$,',...
    ' System ${}^CD_{-\infty}^\alpha x = (a+b\sin(\omega  t))x$',...
    ['Parameter: $a = ', num2str(a), '$, $b = ', num2str(b),...
    '$, $\alpha = ', num2str(alpha), '$, $\omega = ', num2str(omega),...
    '$, $N = ', num2str(N), '$'], ...
    ['Matrix size $2N+1 = ',num2str(2*N+1),'$, number of solutions $\lambda$ found = ', num2str(length(lambdas))]})
hold on
if ~isempty(lambdas)
    help1 = zeros(length(lambdas),1);
    for kk = length(lambdas)
        help1(kk) = log(abs(det(Mat_H_N_alpha(lambdas(kk)))));
    end
    handle_scatter = scatter3(real(lambdas), imag(lambdas), help1, 'ro', 'filled','DisplayName','$\lambda$ found by grid search');
    hLegend =legend(handle_scatter);
end

% Adjust the legend's position within the axes
newPosition = [0.73, 0.3, 0.1, 0.1]; % [x, y, width, height] in normalized units
set(hLegend, 'Position', newPosition);

% % save figure
% savePlot(fig_surf,sprintf('Det_H_N_a_%.2f_b_%.2f_alpha_%.2f_omega_%.2f_N_%d', ...
%        a, b, alpha, omega, N));

%% plot only eigenvalues in complex plane
if ~isempty(lambdas)
    fig_ComplexPlane = figure;
    handle_scatter = scatter(real(lambdas), imag(lambdas), 'ro', 'filled','DisplayName','$\lambda$ found by grid search');
    xlabel('real($\lambda$)')
    ylabel('imag($\lambda$)')

    % Plot the search circle
    theta = linspace(0, 2*pi, 1000);
    x_circle = real(gam) + rad * cos(theta);
    y_circle = imag(gam) + rad * sin(theta);
    hold on;
    plot(x_circle, y_circle, 'b--', 'LineWidth', 1.5, 'DisplayName', ['region with ',num2str(nbr_evs),' estimated \lambda']);
    hold off;

    title({'Solutions of the Fractional Hill Matrix $\det H_N^{(\alpha)}(\lambda)=0$,',...
        ' System ${}^CD_{-\infty}^\alpha x = (a+b\sin(\omega  t))x$',...
        ['Parameter: $a = ', num2str(a), '$, $b = ', num2str(b),...
        '$, $\alpha = ', num2str(alpha), '$, $\omega = ', num2str(omega),...
        '$, $N = ', num2str(N), '$'], ...
        ['Matrix size $2N+1 = $',num2str(2*N+1),', number of solutions $\lambda$ found = ', num2str(length(lambdas))]})

    % Optionally add a legend
    legend('Location', 'best');
end

%% plot magnitudes of vector entries
figure
hold on
for kk = 1:length(lambdas)
p_vec = null(Mat_H_N_alpha(lambdas(kk)));
    p_vec = p_vec/norm(p_vec); 
plot(abs(p_vec))
end
xlabel('entry nr')
ylabel('magnitude of entry')

%% construct the identified time-signal
if ~isempty(lambdas)
    [~, index] = min(abs(imag(lambdas)));
    lambda_sol = lambdas(index);
    
    % DEBUGGING _______________________________________________________________

    % p_vec = X(:,index);
    
    % DEBUGGING END ___________________________________________________________

    p_vec = null(Mat_H_N_alpha(lambda_sol));
    p_vec = p_vec/norm(p_vec); % normalize vector
    omega_vec = (-N:N)'*omega;
    
    soln = @(t) exp(lambda_sol*t).*(p_vec.'*exp(1i*omega_vec*t));
    frac_der_soln = @(t) p_vec.'*((lambda_sol+1i*omega_vec).^alpha.*exp((lambda_sol+1i*omega_vec)*t));
    
    t_vec = -10:0.01:5;

    % soln_compare = zeros(size(t_vec));
    % for kk = 1:length(t_vec)
    %     soln_compare(kk)=soln(t_vec(kk));
    % end

    soln_vec = soln(t_vec);
    frac_der_soln_vec = frac_der_soln(t_vec);
    % DEBUGGING: finding derivative numerically
    t_vec_shifted = t_vec - t_vec(1);
    % frac_der_num_soln_vec = numerical_caputo_derivative(soln_vec, t_vec_shifted, alpha);

    % figure
    % hold on
    % plot(t_vec,real(soln_vec),'DisplayName','$\mathrm{Re}(y(t))$');
    % plot(t_vec,imag(soln_vec),'DisplayName','$\mathrm{Im}(y(t))$');
    % % plot(t_vec,real(soln_compare),'DisplayName','$\mathrm{Re}(y(t))$ compare');
    % % plot(t_vec,imag(soln_compare),'DisplayName','$\mathrm{Im}(y(t))$ compare');   
    % legend
    % figure
    % hold on
    % plot(t_vec,real(frac_der_soln_vec),'DisplayName','$\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % plot(t_vec,imag(frac_der_soln_vec),'DisplayName','$\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    
    syst_mat_fcn = @(t)a +b*sin(t);
    syst_mat_vec = syst_mat_fcn(t_vec); % TBD !
    RHS = syst_mat_vec.*soln_vec;
    % plot(t_vec,real(RHS),'--','DisplayName','$\mathrm{Re}(J(t)y(t))$');
    % plot(t_vec,imag(RHS),'--','DisplayName','$\mathrm{Im}(J(t)y(t))$');
    % legend
    figure
    hold on
    res = frac_der_soln_vec-RHS;
    % plot(t_vec,real(res),'DisplayName','Real res')
    % plot(t_vec,imag(res),'DisplayName','Imag res')
    semilogy(t_vec,abs(res),'DisplayName','residual $|D^\alpha y_{Hill} -J(t)y_{Hill}|$ norm');
    title('residual')
    legend
set(gca, 'YScale', 'log');

    % figure
    % hold on
    % plot(t_vec,abs(soln_vec),'DisplayName','$|(y(t))|$');
    % plot(t_vec,angle(soln_vec),'DisplayName','$\mathrm{arg}(y(t))$');
    % legend
end
%% time integration verification using diffusive repr. time stepping
if ~isempty(lambdas)
    % % assume above section ran, i.e. we have available:
    % % lambda_sol
    % % p_vec
    % % and so on...
    % 
    % % set up the forcing term for infinite memory
    % der_soln = @(t) p_vec.'*((lambda_sol+1i*omega_vec).*exp((lambda_sol+1i*omega_vec)*t));
    % t_L = -Inf; % lower bound for integral
    % forcing = @(t) integral(@(s) (t-s).^(-alpha).*der_soln(s)/gamma(1-alpha),t_L,0);
    % 
    % rhs = @(t,y) syst_mat_fcn(t).*y -forcing(t);
    % rhs_d = @(t,y) syst_mat_fcn(t);
    % 
    % addpath C:\Users\Paul\Desktop\CaputoSolvers
    % Tspan = [0,5];
    % h=0.01;
    % [t_sim, y_sim] = fde12(alpha,rhs,Tspan(1),Tspan(2),soln(0),h);
    % % [t_sim, y_sim]= G1e(alpha,rhs,Tspan,soln(0),h);
    % % [t_sim, y_sim]= G1i(alpha,rhs,Tspan,soln(0),h);
    % % [t_sim, y_sim]= DiffusiveCaputoSolver(alpha, rhs, rhs_d, Tspan, soln(0), 10^4);
    % 
    % soln_Hill = soln(t_sim);
    % 
    % figure
    % hold on
    % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} ~y_{sim}(t)$');
    % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} ~y_{sim}(t)$');
    % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % 
    % plot(t_sim,real(soln_Hill),'--','DisplayName','$\mathrm{Re} ~y_{Hill}(t)$');
    % plot(t_sim,imag(soln_Hill),'--','DisplayName','$\mathrm{Im} ~y_{Hill}(t)$');
    % title('time stepping solution vs Hill')
    % xlabel('$t$')
    % legend
    % set(gca, 'FontSize', 16);
    % 
    % figure
    % diff_sim = y_sim-soln_Hill;
    % rel_diff = abs(diff_sim)./abs(y_sim);
    % % plot(t_sim,abs(real(diff_sim))./abs(y_sim),'DisplayName','Real diff relative')
    % % plot(t_sim,abs(imag(diff_sim))./abs(y_sim),'DisplayName','Imag diff relative')
    % semilogy(t_sim(2:end),rel_diff(2:end),'DisplayName','relative error $|y_{sim}-y_{Hill}|/|y_{sim}|$');
    % set(gca, 'FontSize', 16);
    % xlabel('$t$')
    % legend('Location','southeast')
    % 
    % % DEBUGGING _________________________________________________________
    % % lambda_test = 1+1i;
    % % x_anal = @(t) exp(lambda_test*t);
    % % 
    % % % set up the forcing term for infinite memory
    % % der_x = @(t) x_anal(t)*lambda_test;
    % % t_L = -10; % lower bound for integral
    % % forcing = @(t) integral(@(s) (t-s).^(-alpha).*der_x(s)/gamma(1-alpha),t_L,0);
    % % 
    % % rhs_test = @(t,y) lambda_test^alpha*y -forcing(t);
    % % 
    % % [t_sim, y_sim]= fde12(alpha,rhs_test,0,1,x_anal(0),.001);
    % % 
    % % soln_anal = x_anal(t_sim);
    % % 
    % % figure
    % % hold on
    % % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} y_{sim}(t)$');
    % % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} y_{sim}(t)$');
    % % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % % 
    % % plot(t_sim,real(soln_anal),'--','DisplayName','$\mathrm{Re} y_{anal}(t)$');
    % % plot(t_sim,imag(soln_anal),'--','DisplayName','$\mathrm{Im} y_{anal}(t)$');
    % % title('debugging fde12 exp solution')
    % % legend
    % % 
    % % % MORE DEBUGGING
    % % lambda_test = 1+1i;
    % % z_anal = @(t) mlf(alpha,1,lambda_test*t.^alpha);
    % % 
    % % rhs_test = @(t,y) lambda_test*y;
    % % 
    % % [t_sim, y_sim]= fde12(alpha,rhs_test,0,1,z_anal(0),.001);
    % % 
    % % soln_anal = z_anal(t_sim);
    % % 
    % % figure
    % % hold on
    % % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} y_{sim}(t)$');
    % % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} y_{sim}(t)$');
    % % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % % 
    % % plot(t_sim,real(soln_anal),'--','DisplayName','$\mathrm{Re} y_{anal}(t)$');
    % % plot(t_sim,imag(soln_anal),'--','DisplayName','$\mathrm{Im} y_{anal}(t)$');
    % % title('debugging fde12 mlf solution')
    % % legend
    % % 
    % % % EVEN MORE DEBUGGING
    % % lambda_test = 1+1i;
    % % z_anal = @(t) exp(lambda_test*t);
    % % 
    % % rhs_test = @(t,y) lambda_test*y;
    % % 
    % % [t_sim, y_sim]= ode45(rhs_test,[0,1],z_anal(0));
    % % 
    % % soln_anal = z_anal(t_sim);
    % % 
    % % figure
    % % hold on
    % % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} y_{sim}(t)$');
    % % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} y_{sim}(t)$');
    % % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % % 
    % % plot(t_sim,real(soln_anal),'--','DisplayName','$\mathrm{Re} y_{anal}(t)$');
    % % plot(t_sim,imag(soln_anal),'--','DisplayName','$\mathrm{Im} y_{anal}(t)$');
    % % title('debugging ode45 exp solution')
    % % legend
    % 
    % % DEBUGGING END______________________________________________________
end

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

function df = finiteDiff(f,z,h)
    df = (f(z+h)-f(z-h))/(2*h);
end

function nbr_evs = EstNumEV(f, df, gam, rad, nc)
% EstNumEV estimates the number of roots (eigenvalues) of a complex-valued 
% function f inside a circular contour in the complex plane using the 
% argument principle and the trapezoidal rule.
%
% Inputs:
%   f   - Function handle for the complex-valued function f(z)
%   df  - Function handle for the derivative of f(z), i.e., f'(z)
%   gam - Complex number specifying the center of the circular contour
%   rad - Positive real number specifying the radius of the contour
%   nc  - Positive integer specifying the number of quadrature points 
%         (used for trapezoidal rule integration)
%
% Output:
%   nbr_evs - Estimated number of roots of f inside the specified circle

% Generate unit roots and quadrature points on the contour
w = exp(2i * pi * (1:nc) / nc); 
z = gam + rad * w; 

% Approximate the contour integral using the trapezoidal rule
nr = 0;
for kk = 1:nc
    nr = nr + rad * w(kk) * df(z(kk)) / f(z(kk));
end

% Average the result to estimate the number of eigenvalues
nbr_evs = real(nr / nc); 
end
