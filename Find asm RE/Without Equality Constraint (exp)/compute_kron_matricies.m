% Computes the M matrix such that psi' * Q * psi < 0

function [A1,A2,W,E,II] = compute_kron_matricies(N,M)

% Find differentiation matrices
[D,y] = cheb(N); 

% Square matrix and remove the first row and column due to BC
D2 = D^2; 
D2 = D2(2:N,2:N);

% Creates the diagonal term, ensuring to factor in the BC
S = diag([0; 1 ./(1 - y(2:N).^2); 0]);

% Equation for the fourth derivative matrix
D4 = (diag(1 - y.^2)*D^4 - 8*diag(y)*D^3 - 12*D^2)*S;
D4 = D4(2:N,2:N);

% Sets up the eigenvalue problem by creating two seperate matricies
I = eye(N - 1);

% Creates the elements needed for the Clensure Curtis integration
[~,w] = clencurt(N);
WW = diag(w(2:N));

% Create a 3D tensor to store matricies for a specific n
AA1 = zeros(N-1,N-1,2*M+1);
AA2 = zeros(N-1,N-1,2*M+1);
EE = zeros(N-1,N-1,2*M+1);

% Large matrix to store the matricies along its diagonal
A1 = zeros((N-1)*(2*M+1), (N-1)*(2*M+1));
A2 = zeros((N-1)*(2*M+1), (N-1)*(2*M+1));
E = zeros((N-1)*(2*M+1), (N-1)*(2*M+1));

% Create WW matrix 
W = kron(eye(2*M+1),WW);
II = eye(size(W,1));

% Stores all matricies in a 3D tensor
for n = -M:1:M
    
    % Counter used to store values
    j = n + M + 1;
        
    % Creates the A1 and A2 operator matricies
    AA1(:,:,j) = (-n^4*I + 2*n^2*D2 - D4);   
    AA2(:,:,j) = (2i*I*n + 1i*diag(1 - y(2:N).^2)*(D2*n - I*n^3));
    EE(:,:,j) = (n^2 - D2);
    
    A1(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = AA1(:,:,j);    
    A2(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = AA2(:,:,j); 
    E(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = EE(:,:,j);
            
end
 
% Sets matrix as sparse to speed up operations
% A1 = sparse(A1);
% A2 = sparse(A2);
% W = sparse(W);

end