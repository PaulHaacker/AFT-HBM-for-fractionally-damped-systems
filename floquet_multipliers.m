function lambda = floquet_multipliers(sys_jac, x0p, T, odeopts)
% calculate the Floquet multipliers of a periodic solution
% x0p is a point on the periodic solution, T is the period time

[t,xp,Phi_T] = simulate_x_Phi(sys_jac,[ 0 T],x0p,odeopts);
lambda = eig(Phi_T);

figure(1)
hold on
% plot the Floquet mulipliers lambda in the complex plane
plot(real(lambda),imag(lambda),'xb')

% plot the unit circle
theta = linspace(0,2*pi,100); %to draw the unit circle
plot(real(lambda),imag(lambda),'xb',sin(theta),cos(theta),'r')

axis square
title('Task 2.2: Floquet multipliers')
xlabel('real part')
ylabel('imaginary part')
hold off

end
