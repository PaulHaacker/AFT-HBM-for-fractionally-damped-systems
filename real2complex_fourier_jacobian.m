function J_cell = real2complex_fourier_jacobian(J_realFourierStacked)
% Converts the stacked real Fourier coefficient matrices of J(t) to the
% complex cell array format required by giveFracHill.
%
% Input
%   J_realFourierStacked : n*(2N+1) x n, from jacobian_fourier_coeffs
%
% Output
%   J_cell : cell array of 4N+1 complex n x n matrices
%            ordered as {J_{-2N}, ..., J_0, ..., J_{+2N}}
%            where J_k = 0 for N < |k| <= 2N

n = size(J_realFourierStacked, 2);
N = (size(J_realFourierStacked, 1)/n - 1) / 2;

J_F0 = J_realFourierStacked(1:n, :);
J_FC = permute(reshape(J_realFourierStacked(n+1       : n*(N+1), :), n, N, n), [1,3,2]);
J_FS = permute(reshape(J_realFourierStacked(n*(N+1)+1 : end,     :), n, N, n), [1,3,2]);

% complex coefficients: J_k = (J^Ck - i*J^Sk)/2, J_{-k} = conj(J_k)
J_pos = (J_FC - 1i*J_FS) / 2;   % J_pos(:,:,k) = J_k,   k = 1..N
J_neg = (J_FC + 1i*J_FS) / 2;   % J_neg(:,:,k) = J_{-k}, k = 1..N

zero_block = zeros(n);
J_cell = repmat({zero_block}, 1, 4*N+1);  % initialise with zeros

for k = -N:N
    idx = k + 2*N + 1;  % 1-based index into cell array
    if k == 0
        J_cell{idx} = J_F0;
    elseif k > 0
        J_cell{idx} = J_pos(:,:,k);
    else
        J_cell{idx} = J_neg(:,:,-k);
    end
end
