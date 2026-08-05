% testing the normalzeMatrix function
sysSize = 2;
MatHandle = @(lambda)kron(reshape(1:9,[3 3]),ones(2,2));
alpha = 1;

[norm_MatHandle] = normalizeMatrix(MatHandle,sysSize,alpha);