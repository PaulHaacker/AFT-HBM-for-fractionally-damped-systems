function x_end = simulate_and_plot(sys, x0, t0, tend, odeopts, task_description,fignum,tau)
% simulate the trajectory of sys starting at (t0, x0) until tend. 
% Open a figure, plot the time history (t,x) in the upper subfigure.  
% Plot the phase plane (x1,x2) in the lower subfigure. 
% Set the title of the figure to task_description. 
% Optionally, use tau to scale the time-axis with t_scale in the time-history. 
% Return x(tend).
%
% Input
%   sys                 : handle to a system, specifiying the right-hand side f(t,x)
%   x0                  : initial condition, state at t0
%   t0                  : starting time-instant
%   tend                : end time
%   odeopts             : options for ode45
%   task_description    : string which can be used to set as figure title
%   fignum              : figure number
%   tau                 : period time of the periodic solution
%
% Output
%   x_end               : state at end time tend

if nargin<7, figure, else figure(fignum), end
if nargin<8, 
    t_label = 't';
    t_scale = 1;
else 
    t_label = 't/tau';
    t_scale = tau;
end

[t,x] = ode45(sys,[t0 tend],x0,odeopts);
subplot(2,1,1);
plot(t/t_scale,x)
title(task_description)
xlabel(t_label)
ylabel('x_1, x_2')
legend('x_1','x_2')

subplot(2,1,2);
plot(x(:,1),x(:,2))
xlabel('x_1')
ylabel('x_2')


end
