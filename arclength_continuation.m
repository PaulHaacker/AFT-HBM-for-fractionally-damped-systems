function [mu,x] = arclength_continuation(fun,mu_start,mu_end,x0,tol,maxiter,sigma_start,sigma_max,plot_fun)
% Pseudo-arclength continuation method, see Algorithm 3
%
% Input 
%  fun         : a function defining f(x,mu), the Jacobian df/dx and df/dmu
%  mu_start    : first value of mu
%  mu_end      : last value of mu
%  x0          : initial guess for x which solves f(x,mu_start) = 0
%  tol         : tolerance
%  maxiter     : maximum number of iterations of the Newton method
%  sigma_start : starting value for the stepsize sigma
%  sigma_max   : maximal value for the stepsize sigma
%  plot_fun    : a plotting function
%
% Output
%  mu          : vector with mu values
%  x           : vector with x values

% Remco Leine, INM, University of Stuttgart, 2023

plot_on = nargin>8;

% start up
j = 1;
mu(j) = mu_start;
sigma = sigma_start;
[x(:,j),F,dFdx,conv] = newton(@(x) fun(x,mu_start),x0,tol,maxiter); % find initial point
if ~conv, disp('bad initial point'), end

% find first tangent vector
[F,dFdx,dFdmu] = fun(x(:,j),mu(j));
v_mu = sign(mu_end-mu_start);
v_x = -dFdx\dFdmu*v_mu;
v = [v_x;v_mu];
v = v/norm(v);
v_x = v(1:length(x0));
v_mu = v(end);

if plot_on
    figure, hold on
    xlabel('mu')
    ylabel('x')
end

while (mu(j) <= max(mu_start,mu_end)) && (mu(j)>=min(mu_start,mu_end))

    % predictor
    x_p = x(:,j) + sigma*v_x;
    mu_p = mu(j) + sigma*v_mu;
    if plot_on, plot_fun(mu(j),x(:,j)); drawnow; end

    % corrector
    conv = 0;
    m = 0;
    x_c = x_p;
    mu_c = mu_p;
    while ~conv && m<maxiter
        m = m+1;
        [F,dFdx,dFdmu] = fun(x_c,mu_c);
        H = [dFdx dFdmu;
             v_x'  v_mu];        
        h = [-F;0];
        Dz = H\h;
        Dx = Dz(1:end-1);
        Dmu = Dz(end);
        x_c = x_c + Dx;
        mu_c= mu_c + Dmu;
        if mu_c<=0 || ~isreal(x_c) || ~isreal(mu_c)
            conv = false;   % step overshot into an unphysical (negative/complex) mu
            break
        end
        conv = norm(F)<tol;
    end
    if conv
        j = j+1;
        x(:,j) = x_c;                   % accept new point
        mu(j) = mu_c;
        v = [dFdx dFdmu;
               v']\[zeros(length(x0),1);1];
        v = v/norm(v);                  % compute new tangent
        v_x = v(1:length(x0));
        v_mu = v(end);
        sigma = 1.2*sigma;              % increase stepsize
        if abs(sigma) > sigma_max, sigma = sign(sigma)*sigma_max; end % set maximal stepsize
    else
        sigma = sigma/2;                % reduce stepsize
    end
end
if plot_on
    hold off
end    
