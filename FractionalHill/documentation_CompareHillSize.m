%% Run over array of N and collect results

% % sys parameters
% a = -1;
% b = 2.2;
% alpha = .5;
% omega = 1;
% n = 1;
% 
% % Ns = [5, 10, 15, 20, 25, 30]; % array of N values
% Ns = 5:5:60; % array of N values
% lambdas_all = cell(length(Ns), 1);
% rel_diffs_all = cell(length(Ns), 1);
% t_sim_all = [];
% 
% figure_rel_diff = figure; hold on;
% legend_entries = {};
% 
% % Initialize progress bar
% h_wait = waitbar(0, 'Running simulations...');
% 
% for idx = 1:length(Ns)
%     waitbar(idx/length(Ns), h_wait, sprintf('Running N = %d...', Ns(idx)));
% 
%     N = Ns(idx);
% 
%     % Construct Fourier coefficient vector
%     J = repmat({zeros(n,n)},1,4*N+1);
%     J{2*N} = 1i*b/2; % J_{-1}
%     J{2*N+1} = a;    % J_0
%     J{2*N+2} = -1i*b/2; % J_1
% 
%     % Construct Hill matrix and compute zeros
%     Mat_H_N_alpha = giveFracHill(omega, alpha, J);
%     Mat_norm = normalizeMatrix(Mat_H_N_alpha,length(J{1}),alpha);
%     lambdas = FractionalHillZeros(Mat_norm, x_interval, y_interval, num_points);
% 
%     % Save eigenvalues
%     lambdas_all{idx} = lambdas;
% 
%     if isempty(lambdas), continue; end
% 
%     % Pick dominant lambda and corresponding vector
%     [~, index] = min(abs(imag(lambdas)));
%     lambda_sol = lambdas(index);
%     p_vec = null(Mat_H_N_alpha(lambda_sol));
%     p_vec = p_vec/norm(p_vec);
%     omega_vec = (-N:N)'*omega;
% 
%     soln = @(t) exp(lambda_sol*t).*(p_vec.'*exp(1i*omega_vec*t));
%     der_soln = @(t) p_vec.'*((lambda_sol+1i*omega_vec).*exp((lambda_sol+1i*omega_vec)*t));
%     forcing = @(t) integral(@(s) (t-s).^(-alpha).*der_soln(s)/gamma(1-alpha),-Inf,0);
%     rhs = @(t,y) (a + b*sin(t)).*y - forcing(t);
%     rhs_d = @(t,y) (a + b*sin(t));
% 
%     Tspan = [0,1];
%     h = 0.001;
%     [t_sim, y_sim] = fde12(alpha,rhs,Tspan(1),Tspan(2),soln(0),h);
%     soln_Hill = soln(t_sim);
% 
%     % Store rel_diff
%     diff_sim = y_sim - soln_Hill;
%     rel_diff = abs(diff_sim)./abs(y_sim);
%     rel_diffs_all{idx} = rel_diff;
%     if isempty(t_sim_all)
%         t_sim_all = t_sim;
%     end
% 
%     % Plot on single figure
%     semilogy(t_sim(2:end), rel_diff(2:end), 'DisplayName', sprintf('N = %d', N));
%     legend_entries{end+1} = sprintf('N = %d', N);
% end
% 
% % Close progress bar
% close(h_wait);
% 
% xlabel('$t$', 'Interpreter', 'latex')
% ylabel('relative error', 'Interpreter', 'latex')
% title('Relative difference vs. time for different N', 'Interpreter', 'latex')
% legend(legend_entries, 'Location', 'southeast')
% set(gca, 'FontSize', 14);

%% Run over array of step sizes h and collect results (with ETA and scaling-aware progress)
% System parameters
a = -1;
b = 2.2;
alpha = 0.5;
omega = 1;
n = 1;
N = 30; % fixed N
% hs = [0.01, 0.005, 0.001, 0.0005,0.0001]; % step sizes
hs = logspace(-2,-5,7); % step sizes
rel_diffs_all = cell(length(hs), 1);
t_sim_all = [];
figure_rel_diff = figure; hold on;
legend_entries = {};
% Construct Fourier coefficient vector
J = repmat({zeros(n,n)}, 1, 4*N+1);
J{2*N} = 1i*b/2;     % J_{-1}
J{2*N+1} = a;        % J_0
J{2*N+2} = -1i*b/2;  % J_1
% Hill matrix and eigenvalue analysis
Mat_H_N_alpha = giveFracHill(omega, alpha, J);
Mat_norm = normalizeMatrix(Mat_H_N_alpha, length(J{1}), alpha);
lambdas = FractionalHillZeros(Mat_norm, x_interval, y_interval, num_points);
if isempty(lambdas)
error('No eigenvalues found.');
end
% Use dominant eigenvalue
[~, index] = min(abs(imag(lambdas)));
lambda_sol = lambdas(index);
p_vec = null(Mat_H_N_alpha(lambda_sol));
p_vec = p_vec / norm(p_vec);
omega_vec = (-N:N)' * omega;
soln = @(t) exp(lambda_sol*t) .* (p_vec.' * exp(1i * omega_vec * t));
der_soln = @(t) p_vec.' * ((lambda_sol + 1i * omega_vec) .* exp((lambda_sol + 1i * omega_vec) * t));
forcing = @(t) integral(@(s) (t - s).^(-alpha) .* der_soln(s) / gamma(1 - alpha), -Inf, 0);
rhs = @(t,y) (a + b * sin(t)) .* y - forcing(t);
% Estimate efforts for progress bar
T = 1;
efforts = (1 ./ hs) .* log(1 ./ hs);
efforts = efforts / sum(efforts); % Normalize
cumulative_effort = [0, cumsum(efforts)];
% Initialize waitbar
h_wait = waitbar(0, 'Running simulations...');
% Time the first run to estimate total time
start_time = tic;
h = hs(1);
[t_sim, y_sim] = fde12(alpha, rhs, 0, T, soln(0), h);
elapsed_first = toc(start_time);
estimated_total_time = elapsed_first / efforts(1);
% Store results
soln_Hill = soln(t_sim);
diff_sim = y_sim - soln_Hill;
rel_diff = abs(diff_sim) ./ abs(y_sim);
rel_diffs_all{1} = rel_diff;
t_sim_all = t_sim;
semilogy(t_sim(2:end), rel_diff(2:end), 'DisplayName', sprintf('h = %.4f', h));
legend_entries{end+1} = sprintf('h = %.4f', h);
% Loop over remaining h values
for idx = 2:length(hs)
h = hs(idx);
sim_start = tic;
[t_sim, y_sim] = fde12(alpha, rhs, 0, T, soln(0), h);
soln_Hill = soln(t_sim);
diff_sim = y_sim - soln_Hill;
rel_diff = abs(diff_sim) ./ abs(y_sim);
rel_diffs_all{idx} = rel_diff;
semilogy(t_sim(2:end), rel_diff(2:end), 'DisplayName', sprintf('h = %.4f', h));
legend_entries{end+1} = sprintf('h = %.4f', h);
% Update progress and ETA
progress = cumulative_effort(idx+1);
elapsed_so_far = toc(start_time);
est_remaining = estimated_total_time - elapsed_so_far;
waitbar(progress, h_wait, ...
sprintf('Running h = %.4f... ETA: %.1f sec', h, est_remaining));
end
close(h_wait);
legend show;
xlabel('Time');
ylabel('Relative Difference');
title('Relative Difference for Various Step Sizes h, N = 25');
grid on;
set(gca, 'YScale', 'log');

%% debugging plot eigenvalues

% % Plot stored lambdas in the complex plane for each N
% figure;
% num_N = length(Ns);
% num_cols = ceil(sqrt(num_N));     % number of columns in subplot grid
% num_rows = ceil(num_N / num_cols); % number of rows in subplot grid
% for idx = 1:num_N
% lambdas = lambdas_all{idx};
% subplot(num_rows, num_cols, idx);
% if ~isempty(lambdas)
% plot(real(lambdas), imag(lambdas), 'bo', 'MarkerFaceColor', 'b');
% xline(0, '--k');
% yline(0, '--k');
% grid on;
% title(sprintf('N = %d', Ns(idx)));
% xlabel('Re(\lambda)');
% ylabel('Im(\lambda)');
% else
% title(sprintf('N = %d (no data)', Ns(idx)));
% axis off;
% end
% end
% sgtitle('Eigenvalues in the Complex Plane for Different N');