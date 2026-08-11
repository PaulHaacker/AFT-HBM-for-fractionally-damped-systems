function Mat_H_N_alpha = giveHill(omega, alpha, J, C0)
    % function handing back the Hill matrix for a system with fractional damping % given as \dot y = J(t)y + C0 D^{\alpha}y, where J(t) is a periodic matrix with period T=2*pi/omega and D^{\alpha} is the fractional derivative of order alpha.
    % takes input
    % omega -   angular frequency of the periodic Matrix J, i.e. period
    %           T=2*pi/omega
    % alpha -   fractional order of the fractional damping term
    % J -       Fourier coefficients of J(t) \in R^{n \times n}, need to be
    %           a cell array of 4N+1 entries, each \in C^{n \times n},
    %           ordered as J = {J_{-2N},...,J_0,...,J_{+2N}}
    % C0 -      constant coefficient of the damping term C0 \in R^{n \times  
    % %           n} (i.e. it is not time-dependent)
    test = size(cell2mat(J));
    n = test(1);           % n - system size
    N = (test(2)/n-1)/4;   % N - number of frequencies considered in truncation
    if N~=int8(N)
        error('wrong dimension of J')
    end
    mat_indx = toeplitz(2*N+1:1:4*N+1,2*N+1:-1:1);
    help1 = cell2mat(J(mat_indx));
    Mat_H_N_alpha = @(lambda) help1+kron(diag((lambda+1i*omega*(-N:N)).^alpha),C0)...
        +kron(diag(1i*omega*(N:-1:-N)),eye(n)) -lambda*eye(n*(2*N+1));
end