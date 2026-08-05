% testing giveFracHill
alpha = 1;
omega = 1;
J = {-2*ones(2),-1*ones(2),zeros(2),1*ones(2),2*ones(2)};

Mat_H_N_alpha = giveFracHill(omega, alpha, J);

Mat_H_N_alpha(0)