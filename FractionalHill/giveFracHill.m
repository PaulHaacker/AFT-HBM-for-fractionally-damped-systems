function Mat_H_N_alpha = giveFracHill(omega, alpha, J)
    % takes input
    % omega -   angular frequency of the periodic Matrix J, i.e. period
    %           T=2*pi/omega
    % alpha -   fractional order of the diff eqn
    % J -       Fourier coefficients of J(t) \in R^{n \times n}, need to be
    %           a cell array of 4N+1 entries, each \in C^{n \times n},
    %           ordered as J = {J_{-2N},...,J_0,...,J_{+2N}}
    test = size(cell2mat(J));
    n = test(1);           % n - system size
    N = (test(2)/n-1)/4;   % N - number of frequencies considered in truncation
    if N~=int8(N)
        error('wrong dimension of J')
    end
    mat_indx = toeplitz(2*N+1:1:4*N+1,2*N+1:-1:1);
    help1 = cell2mat(J(mat_indx));
    Mat_H_N_alpha = @(lambda) help1-kron(diag((lambda+1i*omega*(-N:N)).^alpha),eye(n));
end