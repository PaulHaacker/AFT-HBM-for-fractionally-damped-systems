function zdot = sys_ext(t,z,sys_jac)
% system in extended form z' = f_ext(t,z)
% with the extended state z = [x1 .. xn Phi11 ... Phinn]^T
%
% Input
%  t       : time instant
%  z       : extended state
%  sys_jac : system with right-hand side and Jacobian
%
% Output
%  zdot : time-derivative of z

% Remco Leine, INM, University of Stuttgart, 2023
N = length(z);
n = (-1+sqrt(1+4*N))/2;
x = z(1:n);
Phi = vec2mat(z(n+1:end));
[f,dfdx] = feval(sys_jac,t,x);
Phidot = dfdx*Phi;
zdot = [f;mat2vec(Phidot)];
