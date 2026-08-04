function plot_fun_HBM(omega,X,n,N,L)
[t,x] = fourier_coeff_to_time_series(X,omega,n,N,L);
A = max(x(:,1));
plot(omega,A,'or-')
xlabel('omega')
ylabel('max q(t)')