%% script applying the fractional hill method to a minimal example system
% considering the system D^\alpha x = A(t)x(t), fractional Mathieu

clear
close all

%% parameter stash

alpha = 0.7;
a = -1;
b = 5;
omega = 1;
n = 2;
syst_mat_fcn = @(t) [0,1;a+b*sin(omega*t),0];
fdefun = @(t,x) [0,1;a+b*sin(omega*t),0]*x;

%% calling FractionalHillZeros

N = 20; % (2N+1) is the number of fourier coefficients considered

% % search grid in the complex plane
% x_interval = [-5, 120]; % real interval for grid search
% y_interval = [-20, 20]; % imaginary interval for grid search
% num_points = 300; % Number of grid points to use for the search in each dimension.

x_interval = [-1, .25];[-1, 1]; % real interval for grid search
y_interval = [0,5]; % imaginary interval for grid search
num_points = 20; % Number of grid points to use for the search in each dimension.

% x_interval = [0,.1081]; % real interval for grid search
% y_interval = [0,1]; % imaginary interval for grid search
% num_points = 2; % Number of grid points to use for the search in each dimension.

% % Fourier coefficients are available analytically
% J = cell(1,4*N+1); % create cell with zeros
J = repmat({zeros(n,n)},1,4*N+1);
J{2*N} = [0, 0; 1i*b/2, 0]; % J_{-1}
J{2*N+1} = [0, 1; a, 0]; % J_0
J{2*N+2} = [0, 0; -1i*b/2, 0]; % J_1

tic
Mat_H_N_alpha = giveFracHill(omega, alpha, J);
Mat_norm = normalizeMatrix(Mat_H_N_alpha,length(J{1}),alpha);
% lambdas =[];
lambdas = FractionalHillZeros(Mat_norm, x_interval, y_interval, num_points);
toc

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

%% test plot of the found lambda
% r and imaginary grids
grid_size_real = 100;
grid_size_imag = 100;
real_grid = linspace(x_interval(1), x_interval(2), grid_size_real);
imag_grid = linspace(y_interval(1), y_interval(2), grid_size_imag);
[real_mesh, imag_mesh] = meshgrid(real_grid, imag_grid);

% Compute |det M(lambda)|
absDet_mesh = zeros(size(real_mesh));
for kk = 1:length(real_grid)
    for jj = 1:length(imag_grid)
        % absDet_mesh(jj, kk) = abs(det(Mat_H_N_alpha(real_grid(kk) + 1i * imag_grid(jj))));
        absDet_mesh(jj, kk) = abs(det(Mat_norm(real_grid(kk) + 1i * imag_grid(jj))));
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
    title({'Solutions of the Fractional Hill Matrix $\det H_N^{(\alpha)}(\lambda)=0$,',...
        ' System ${}^CD_{-\infty}^\alpha x = (a+b\sin(\omega  t))x$',...
        ['Parameter: $a = ', num2str(a), '$, $b = ', num2str(b),...
        '$, $\alpha = ', num2str(alpha), '$, $\omega = ', num2str(omega),...
        '$, $N = ', num2str(N), '$'], ...
        ['Matrix size $2N+1 = $',num2str(2*N+1),', number of solutions $\lambda$ found = ', num2str(length(lambdas))]})
    % save figure
    % savePlot(fig_ComplexPlane,sprintf('Det_H_N_a_%.2f_b_%.2f_alpha_%.2f_omega_%.2f_N_%d_ComplexPlane', ...
    %        a, b, alpha, omega, N));
end

lambdas_all = lambdas;
lambdas = lambdas(real(lambdas) > 0);
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
    p_vec = reshape(p_vec,[n,2*N+1]);
    omega_vec = (-N:N)'*omega;

    soln = @(t) exp(lambda_sol*t).*(p_vec*exp(1i*omega_vec*t));
    frac_der_soln = @(t) p_vec*((lambda_sol+1i*omega_vec).^alpha.*exp((lambda_sol+1i*omega_vec)*t));

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

    % syst_mat_fcn = @(t)a +b*sin(t);
    RHS = zeros(size(soln_vec));
    for kk = 1:length(t_vec)
        syst_mat_stage = syst_mat_fcn(t_vec(kk)); % TBD !
        RHS(:,kk) = syst_mat_stage*soln_vec(:,kk);
    end
    % plot(t_vec,real(RHS),'--','DisplayName','$\mathrm{Re}(J(t)y(t))$');
    % plot(t_vec,imag(RHS),'--','DisplayName','$\mathrm{Im}(J(t)y(t))$');
    % legend
    figure
    hold on
    res = frac_der_soln_vec-RHS;
    % plot(t_vec,real(res),'DisplayName','Real res')
    % plot(t_vec,imag(res),'DisplayName','Imag res')
    plot(t_vec,abs(res),'DisplayName','residual $|D^\alpha y_{Hill} -J(t)y_{Hill}|$ norm');
    set(gca,'YScale','log')
    title('residual')
    legend

    % figure
    % hold on
    % plot(t_vec,abs(soln_vec),'DisplayName','$|(y(t))|$');
    % plot(t_vec,angle(soln_vec),'DisplayName','$\mathrm{arg}(y(t))$');
    % legend
end
%% time integration verification using diffusive repr. time stepping
if ~isempty(lambdas)
    % assume above section ran, i.e. we have available:
    % lambda_sol
    % p_vec
    % and so on...

    % set up the forcing term for infinite memory
    der_soln = @(t) p_vec*((lambda_sol+1i*omega_vec).*exp((lambda_sol+1i*omega_vec)*t));
    t_L = -Inf; % lower bound for integral
    forcing = @(t) integral(@(s) (t-s).^(-alpha).*der_soln(s)/gamma(1-alpha),t_L,0,'ArrayValued',true);

    rhs = @(t,y) syst_mat_fcn(t)*y -forcing(t);
    rhs_d = @(t,y) syst_mat_fcn(t);

    addpath C:\Users\Paul\Desktop\CaputoSolvers
    Tspan = [0,1];
    h=0.01;
    [t_sim, y_sim] = fde12(alpha,rhs,Tspan(1),Tspan(2),soln(0),h);
    % [t_sim, y_sim]= G1e(alpha,rhs,Tspan,soln(0),h);
    % [t_sim, y_sim]= G1i(alpha,rhs,Tspan,soln(0),h);
    % [t_sim, y_sim]= DiffusiveCaputoSolver(alpha, rhs, rhs_d, Tspan, soln(0), 10^4);

    soln_Hill = soln(t_sim);

    figure
    hold on
    plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} ~y_{sim}(t)$');
    plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} ~y_{sim}(t)$');
    % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');

    plot(t_sim,real(soln_Hill),'--','DisplayName','$\mathrm{Re} ~y_{Hill}(t)$');
    plot(t_sim,imag(soln_Hill),'--','DisplayName','$\mathrm{Im} ~y_{Hill}(t)$');
    title('time stepping solution vs Hill')
    xlabel('$t$')
    legend
    set(gca, 'FontSize', 16);

    figure
    diff_sim = y_sim-soln_Hill;
    rel_diff = vecnorm(diff_sim)./vecnorm(y_sim);
    % plot(t_sim,abs(real(diff_sim))./abs(y_sim),'DisplayName','Real diff relative')
    % plot(t_sim,abs(imag(diff_sim))./abs(y_sim),'DisplayName','Imag diff relative')
    semilogy(t_sim(2:end),rel_diff(2:end),'DisplayName','relative error $|y_{sim}-y_{Hill}|/|y_{sim}|$');
    set(gca, 'FontSize', 16);
    xlabel('$t$')
    legend('Location','southeast')

    % DEBUGGING _________________________________________________________
    % lambda_test = 1+1i;
    % x_anal = @(t) exp(lambda_test*t);
    % 
    % % set up the forcing term for infinite memory
    % der_x = @(t) x_anal(t)*lambda_test;
    % t_L = -10; % lower bound for integral
    % forcing = @(t) integral(@(s) (t-s).^(-alpha).*der_x(s)/gamma(1-alpha),t_L,0);
    % 
    % rhs_test = @(t,y) lambda_test^alpha*y -forcing(t);
    % 
    % [t_sim, y_sim]= fde12(alpha,rhs_test,0,1,x_anal(0),.001);
    % 
    % soln_anal = x_anal(t_sim);
    % 
    % figure
    % hold on
    % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} y_{sim}(t)$');
    % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} y_{sim}(t)$');
    % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % 
    % plot(t_sim,real(soln_anal),'--','DisplayName','$\mathrm{Re} y_{anal}(t)$');
    % plot(t_sim,imag(soln_anal),'--','DisplayName','$\mathrm{Im} y_{anal}(t)$');
    % title('debugging fde12 exp solution')
    % legend
    % 
    % % MORE DEBUGGING
    % lambda_test = 1+1i;
    % z_anal = @(t) mlf(alpha,1,lambda_test*t.^alpha);
    % 
    % rhs_test = @(t,y) lambda_test*y;
    % 
    % [t_sim, y_sim]= fde12(alpha,rhs_test,0,1,z_anal(0),.001);
    % 
    % soln_anal = z_anal(t_sim);
    % 
    % figure
    % hold on
    % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} y_{sim}(t)$');
    % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} y_{sim}(t)$');
    % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % 
    % plot(t_sim,real(soln_anal),'--','DisplayName','$\mathrm{Re} y_{anal}(t)$');
    % plot(t_sim,imag(soln_anal),'--','DisplayName','$\mathrm{Im} y_{anal}(t)$');
    % title('debugging fde12 mlf solution')
    % legend
    % 
    % % EVEN MORE DEBUGGING
    % lambda_test = 1+1i;
    % z_anal = @(t) exp(lambda_test*t);
    % 
    % rhs_test = @(t,y) lambda_test*y;
    % 
    % [t_sim, y_sim]= ode45(rhs_test,[0,1],z_anal(0));
    % 
    % soln_anal = z_anal(t_sim);
    % 
    % figure
    % hold on
    % plot(t_sim,real(y_sim),'DisplayName','$\mathrm{Re} y_{sim}(t)$');
    % plot(t_sim,imag(y_sim),'DisplayName','$\mathrm{Im} y_{sim}(t)$');
    % % plot(t_vec,real(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Re}\mathrm{D}^\alpha y(t)$');
    % % plot(t_vec,imag(frac_der_num_soln_vec),'DisplayName','numerical $\mathrm{Im}\mathrm{D}^\alpha y(t)$');
    % 
    % plot(t_sim,real(soln_anal),'--','DisplayName','$\mathrm{Re} y_{anal}(t)$');
    % plot(t_sim,imag(soln_anal),'--','DisplayName','$\mathrm{Im} y_{anal}(t)$');
    % title('debugging ode45 exp solution')
    % legend

    % DEBUGGING END______________________________________________________
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