% Comparing the use of sparse and not sparse matricies
clear all

% Set discritisation grid
N = 25;
M = 25;

% Set Reynolds number
Re = 87.7173;
lamda = 1/Re;

% Call matricies
[A1,A2,W] = compute_kron_matricies(N,M);

% Create sparse matricies
AS1 = sparse(A1);
AS2 = sparse(A2);
WS = sparse(W);

% Compute matrix M
M = W*(lamda*A1 + A2) + (lamda*A1 + A2)'*W';
MS = WS*(lamda*AS1 + AS2) + (lamda*AS1 + AS2)'*WS';

% Compute eigenvalues
e = eig(M);
eS = eigs(MS,size(MS,2));

diff_e = e - eS;
max_diff_e = max(diff_e);

