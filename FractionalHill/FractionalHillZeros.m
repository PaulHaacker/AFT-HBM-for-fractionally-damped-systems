function lambdas = FractionalHillZeros(Mat_H_N_alpha, x_interval, y_interval, num_points, varargin)
    % FractionalHillZeros finds zeros of a matrix-valued function in a complex rectangle.
    %
    % Required:
    % Mat_H_N_alpha: Function handle to matrix-valued function.
    % x_interval: [x_min, x_max] real part of search domain.
    % y_interval: [y_min, y_max] imaginary part of search domain.
    % num_points: Number of grid points in each direction.
    %
    % Optional:
    % method: (string) Solver method, default: 'NewtonComplexScalarVectorized'
    
    % Handle optional input
    if nargin >= 5
        method = varargin{1};
    else
        method = 'NewtonComplexScalarVectorized';
    end

 switch lower(method)
    case 'newtoncomplexscalarvectorized'
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-15, 'TolX', 1e-15);
        x_vals = linspace(x_interval(1), x_interval(2), num_points);
        y_vals = linspace(y_interval(1), y_interval(2), num_points);

        lambdas = [];
        zero_fcn_base = @(x) vectorizeComplex(Mat_H_N_alpha, x);
        zero_fcn = zero_fcn_base;
        deflation_terms = {};  % Store previous solutions
        tolerance_duplicates = 1e-4;

        for ii = 1:num_points
            for jj = 1:num_points
                x_0 = [x_vals(ii); y_vals(jj)];
                [x_solved, ~, exitflag] = newton_numericalGradient(zero_fcn, x_0, 1e-13, 100);
                if exitflag > 0
                    lambda = x_solved(1) + 1i * x_solved(2);
                    diffs = abs(lambda - lambdas);
                    if all(diffs > tolerance_duplicates)
                        deflation_terms{end+1} = x_solved;
                        zero_fcn = @(x) zero_fcn_base(x) .* deflation_multiplier(x, deflation_terms);
                        lambdas = [lambdas; lambda];
                    end
                end
            end
        end

    case 'newtoncomplexscalar'
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-15, 'TolX', 1e-15);
        x_vals = linspace(x_interval(1), x_interval(2), num_points);
        y_vals = linspace(y_interval(1), y_interval(2), num_points);

        lambdas = [];
        zero_fcn_base = @(x) det(Mat_H_N_alpha(x));
        zero_fcn = zero_fcn_base;
        deflation_terms = {};  % Store previous solutions
        tolerance_duplicates = 1e-4;

        for ii = 1:num_points
            for jj = 1:num_points
                x_0 = x_vals(ii)+1i* y_vals(jj);
                [lambda, ~, exitflag] = newton_numericalGradient(zero_fcn, x_0, 1e-13, 100);
                if exitflag > 0
                    diffs = abs(lambda - lambdas);
                    if all(diffs > tolerance_duplicates)
                        deflation_terms{end+1} = lambda;
                        zero_fcn = @(x) zero_fcn_base(x) .* deflation_multiplier(x, deflation_terms);
                        lambdas = [lambdas; lambda];
                    end
                end
            end
        end

    otherwise
        error('Unsupported method: %s', method);
end


    % if ~isempty(lambdas)
    %     % Filter duplicates (considering numerical precision)
    %     tolerance_duplicates = 1e-4;
    %     filtered_lambdas = lambdas(1);
    %     % Iterate through the rest of the 'lambdas' vector
    %     for i = 2:length(lambdas)
    %         % Calculate the difference between the current lambda and all the filtered ones
    %         diffs = abs(filtered_lambdas - lambdas(i));
    %         % Check if all differences are greater than the tolerance
    %         if all(diffs > tolerance_duplicates)
    %             % If true, add the current lambda to the filtered set
    %             filtered_lambdas = [filtered_lambdas; lambdas(i)];
    %         end
    %     end
    %     % Output the filtered set
    %     lambdas = filtered_lambdas;
    % end
end

% function lambdas = FractionalHillZeros(Mat_H_N_alpha, x_interval, y_interval, num_points)
%     % Find all zeros of a function in a given complex rectangle interval.
%     %
%     % Mat_H_N_alpha: Function handle to matrix-valued function.
%     % x_interval: [x_min, x_max] for the real part of lambda.
%     % y_interval: [y_min, y_max] for the imaginary part of lambda.
%     % num_points: Number of grid points to use for the search in each dimension.
%     %
%     % Example call:
%     % alpha = 1;
%     % omega = 1;
%     % J = {-2*ones(2),-1*ones(2),zeros(2),1*ones(2),2*ones(2)};
%     % Mat_H_N_alpha = giveFracHill(omega, alpha, J);
%     % x_interval = [-10, 10];
%     % y_interval = [-10, 10];
%     % num_points = 100;
%     % lambdas = FractionalHillZeros(Mat_H_N_alpha, x_interval, y_interval, num_points);
% 
%     options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-15, 'TolX', 1e-15);
%     x_vals = linspace(x_interval(1), x_interval(2), num_points);
%     y_vals = linspace(y_interval(1), y_interval(2), num_points);
% 
%     lambdas = [];
%     % for deflation:
%     zero_fcn_base = @(x) vectorizeComplex(Mat_H_N_alpha, x);
%     zero_fcn = zero_fcn_base;
%     deflation_terms = {};  % Store previous solutions
% 
%     % Filter duplicates (considering numerical precision)
%     tolerance_duplicates = 1e-4;
% 
%     for ii = 1:num_points
%         for jj = 1:num_points
%             x_0 = [x_vals(ii); y_vals(jj)];
%             % [x_solved, ~, exitflag] = fsolve(@(x) vectorizeComplex(Mat_H_N_alpha, x), x_0, options);
%             [x_solved, ~, exitflag] = newton_raphson_numerical( zero_fcn, x_0, 1e-13, 100);
%             if exitflag > 0 % fsolve successfully found a solution
%                 lambda = x_solved(1) + 1i * x_solved(2);
%                 % zero_fcn = @(x) zero_fcn(x)*norm(x-x_solved)^-1;
% 
%                 % check if new lambda has been found
%                 diffs = abs(lambda - lambdas);
%                 % Check if all differences are greater than the tolerance
%                 if all(diffs > tolerance_duplicates)
%                     % If true, add the current lambda to the filtered set
%                     deflation_terms{end+1} = x_solved;  % Store solution
%                     zero_fcn = @(x) zero_fcn_base(x) .* deflation_multiplier(x, deflation_terms);
%                     lambdas = [lambdas; lambda];
%                 end
%             end
%         end
%     end
% 
% end

function out = vectorizeComplex(M, x)
    % Takes x, which is 1x2, representing a complex number lambda with x(1) and
    % x(2) being the real and imaginary part, respectively.
    % Returns the determinant of M(lambda), again vectorized to real and
    % imaginary part.
    M_val = det(M(x(1) + 1i * x(2)));
    out = [real(M_val); imag(M_val)];
    % out = log(abs([real(M_val), imag(M_val)])+1);
    % out = abs(M_val);
end

% function [root, iter, exitflag] = newton_raphson_numerical(f, x0, tol, max_iter)
% % newton_raphson_numerical finds the root of a function using Newton-Raphson
% % with numerical derivative.
% %
% % Inputs:
% %   f         - function handle to the function f(x)
% %   x0        - initial guess
% %   tol       - tolerance for convergence
% %   max_iter  - maximum number of iterations
% %
% % Outputs:
% %   root      - estimated root of f(x) = 0
% %   iter      - number of iterations performed

%     h = 1e-8; % step size for numerical derivative
%     x = x0;
%     iter = 0;

%     % assume quadratic Jacobian
%     n = length(x0);
%     J = zeros(n, n);

%     while iter < max_iter
%         fx = f(x);
  
%         % Estimate Jacobian numerically 
%         for j = 1:n
%             dx = zeros(n,1);
%             dx(j) = h;
%             J(:,j) = (f(x + dx) - f(x - dx)) / (2*h); % (central differences)
%             % J(:,j) = imag(f(x + 1i*dx) ) / (h); % (complex step)
%         end
%         % dfx = (f(x + h) - f(x - h)) / (2 * h); % central difference

%         % if abs(dfx) < eps
%         %     error('Derivative is too small. No convergence.');
%         % end

%         if rcond(J)<1e-10
%             exitflag = 0;
%             root = [];
%             return;
%         end

%         x_new = x - J\fx;
%         % if isnan(x_new(2)) || isnan(x_new(2))
%         %     1;
%         % end
%         if abs(x_new - x) < tol
%             root = x_new;
%             exitflag = 1;
%             return;
%         end

%         x = x_new;
%         iter = iter + 1;
%     end

%     exitflag = 0;
%     root = [];
% end

function mult = deflation_multiplier(x, roots)
    mult = 1;
    for k = 1:length(roots)
        mult = mult * norm(x - roots{k})^-1;
    end
end
