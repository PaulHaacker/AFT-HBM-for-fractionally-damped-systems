function v = mat2vec(M)
% mat2vec transforms a square matrix to a supervector
% M = [v1 v2 v3] -> v = [v1; v2; v3]
n = size(M,1);
v = reshape(M',n^2,1);