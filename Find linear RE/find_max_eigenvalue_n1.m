% Calculates the maximum eigenvalue for a given Reynolds number using the
% Orr-Sommerfeld equation for n = 1

function [max_eig] = find_max_eigenvalue_n1(Re)

% Set discritisation grid to give large accuracy
N = 100;

% Find differentiation matrices
[D,x] = cheb(N);

% Square matrix and remove the first row and column due to BC
D2 = D^2; 
D2 = D2(2:N,2:N);

% Creates the diagonal term, ensuring to factor in the BC
S = diag([0; 1 ./(1 - x(2:N).^2); 0]);

% Equation for the fourth derivative matrix
D4 = (diag(1 - x.^2)*D^4 - 8*diag(x)*D^3 - 12*D^2)*S;
D4 = D4(2:N,2:N);

% Sets up the eigenvalue problem by creating two seperate matricies
I = eye(N - 1);
A = (D4 - 2*D2 + I)/Re - 2i*I - 1i*diag(1 - x(2:N).^2)*(D2 - I);
B = D2 - I;

% Calculates the eigenvalues from the two matrices
ee = eig(A,B);

% Finds the maximum real part of the eigenvalues
max_eig = max(real(ee));

end