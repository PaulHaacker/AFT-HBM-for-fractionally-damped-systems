%% MATLAB Workshop Periodic Solutions
%% Part 3: The harmonic balance method

odeopts = odeset('RelTol',1e-5,'AbsTol',1e-7);

%% Task 3.1.1:
omega = 1;
n = 1;
N = 2;
L = 100;
[V,W,D] = fourier_matrices(omega,n,N,L);
disp('If your implementation is correct, then this should give (almost) zero:')
norm(W'*V - eye(n*(2*N+1)))

%% Task 3.1.2:
omega = 1;
figure(1);
h = gcf; h.Name = 'Task 3.1.2';
fun = @(t) 1  + 2*cos(omega*t) + 3*sin(2*omega*t) + 1/4*cos(5*omega*t);
fun_derivative = @(t) -2*omega*sin(omega*t) + 3*2*omega*cos(2*omega*t) - 5/4*omega*sin(5*omega*t);
t_exact = linspace(0,2*pi/omega,1000);
x_exact = fun(t_exact); 
xdot_exact = fun_derivative(t_exact);
for casenum=1:4
    switch casenum
        case 1
            N = 5; L = 100;
        case 2
            N = 2; L = 100;
        case 3
            N = 5; L = 11;
        case 4
            N = 5; L = 9;
    end  

    [t,x,xdot,x_proj,x_proj_dot] = fourier_test(fun,fun_derivative,omega,1,N,L);
    subplot(2,4,casenum);
    plot(t_exact,x_exact,'b',t,x,'k.',t,x_proj,'rx')
    title(['case ',num2str(casenum),'), N = ',num2str(N),', L = ', num2str(L)])
    xlabel('t')
    legend('x (cont.)','x (sampled)','x_{proj}')

    subplot(2,4,casenum+4);
    plot(t_exact,xdot_exact,'b',t,xdot,'k.',t,x_proj_dot,'rx')
    xlabel('t')
    legend('xdot (cont.)','xdot (sampled)','xdot_{proj}')
end



%% Task 3.2
alpha = 1;
beta = 0.04;
gamma = 1;
delta = 0.1;
omega = 1.3;
tau = 2*pi/omega;
[t,x] = ode45(@(t,x) duffing(t,x,alpha,beta,gamma,delta,omega),[0 100*tau],[0 0],odeopts);
x0p = x(end,:)';

L = 100;
N = 5;
n = 2;
tol = 1e-7;
maxiter = 100;
T = 2*pi/omega;
dt = T/L;
t = (0:dt:T-dt)';
[t_exact,x_exact] = ode45(@(t,x) duffing(t,x,alpha,beta,gamma,delta,omega),t,x0p,odeopts);

X0 = zeros(n*(2*N+1),1);
[X,r,conv,k] = HBM_nonaut(@(t,x) duffing_jac(t,x,alpha,beta,gamma,delta,omega),X0,omega,n,N,L,tol,maxiter);
[t,x] = fourier_coeff_to_time_series(X,omega,n,N,L);

figure(2);
h = gcf; h.Name = 'Task 3.2';
plot(t_exact,x_exact,t,x,'o')
title('HBM solution of the forced Duffing system')
xlabel('t')
ylabel('q, qdot')
legend('q exact','qdot exact','q HBM','qdot HBM')

%% Task 3.3
[mu,x] = arclength_continuation(@(x,mu) truss(x,mu,1,1.2,1,1,0.1),-2,2,[0;0],1e-5,100,1e-2,2e-1);
figure(3);
h = gcf; h.Name = 'Task 3.3';
plot(mu,x(1,:),'o-')
title('Truss system')
xlabel('F')
ylabel('q')

%% Task 3.4
alpha = 1;
beta = 0.04;
gamma = 1;
delta = 0.0;0.1;
omega_start = 2;0.1;
omega_end = 0.1;2;
L = 100;
N = 5;
n = 2;
X0 = zeros(n*(2*N+1),1);

[X0,r,conv,k] = HBM_nonaut(@(t,x) duffing_jac(t,x,alpha,beta,gamma,delta,omega),X0,omega_start,n,N,L,1e-5,100);
% [om,X] = arclength_continuation_HBM(@(t,x,omega) duffing_jac(t,x,alpha,beta,gamma,delta,omega),0.1,2,X0,n,N,L,1e-5,100,1e-2,0.1,@(om,X) plot_fun_HBM(om,X,n,N,L));
[om,X] = arclength_continuation_HBM(@(t,x,omega) duffing_jac(t,x,alpha,beta,gamma,delta,omega),omega_start,omega_end,X0,n,N,L,1e-5,100,1e-2,0.1);

A = zeros(length(om),1);
for i=1:length(om)
  [t,x] = fourier_coeff_to_time_series(X(:,i),om(i),n,N,L);
  A(i) = max(x(:,1));
end
plot(om,A,'o-')
xlabel('$\omega$','Interpreter','latex')
ylabel('$\textrm{max} \, q(t)$','Interpreter','latex')
title('Nonlinear FRF of the Duffing system','Interpreter','latex')
h = gcf; h.Name = 'Task 3.4';