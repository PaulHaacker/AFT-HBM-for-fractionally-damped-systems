function [norm_MatHandle] = normalizeMatrix(MatHandle,sysSize,alpha)
% Mathandle - function handle with scalar, complex input returning complex
% matrix of square size n*(2N+1)
% sysSize - original system size, above denoted n
% alpha - fractional system order

% first extract N
N = 1/2*(length(MatHandle(0))/sysSize -1);

vec1 = (1./(N:-1:1)/sysSize).^alpha;
vec2 = [vec1,1,fliplr(vec1)]';
mat1 = repmat(vec2,1,2*N+1);
mat2 = kron(mat1,ones(sysSize,sysSize));
norm_MatHandle = @(lambda)mat2.*MatHandle(lambda);
end