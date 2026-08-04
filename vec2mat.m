function M = vec2mat(v)
% vec2mat transforms a supervector to a square matrix
% v = [v1; v2; v3]    ->  M = [v1 v2 v3]
% The length of v must be n^2
n = sqrt(length(v));
M = reshape(v(:),n,n)';