function [root, iter, exitflag] = newton_numericalGradient(f, x0, tol, max_iter)
% newton_raphson_numerical finds the root of a function using Newton-Raphson
% with numerical derivative.
%
% Inputs:
%   f         - function handle to the function f(x)
%   x0        - initial guess
%   tol       - tolerance for convergence
%   max_iter  - maximum number of iterations
%
% Outputs:
%   root      - estimated root of f(x) = 0
%   iter      - number of iterations performed

    h = 1e-8; % step size for numerical derivative
    x = x0;
    iter = 0;

    % assume quadratic Jacobian
    n = length(x0);
    J = zeros(n, n);

    while iter < max_iter
        fx = f(x);
  
        % Estimate Jacobian numerically 
        for j = 1:n
            dx = zeros(n,1);
            dx(j) = h;
            J(:,j) = (f(x + dx) - f(x - dx)) / (2*h); % (central differences)
            % J(:,j) = imag(f(x + 1i*dx) ) / (h); % (complex step)
        end
        % dfx = (f(x + h) - f(x - h)) / (2 * h); % central difference

        % if abs(dfx) < eps
        %     error('Derivative is too small. No convergence.');
        % end

        if rcond(J)<1e-10
            exitflag = 0;
            root = [];
            return;
        end

        x_new = x - J\fx;
        % if isnan(x_new(2)) || isnan(x_new(2))
        %     1;
        % end
        if abs(x_new - x) < tol
            root = x_new;
            exitflag = 1;
            return;
        end

        x = x_new;
        iter = iter + 1;
    end

    exitflag = 0;
    root = [];
end