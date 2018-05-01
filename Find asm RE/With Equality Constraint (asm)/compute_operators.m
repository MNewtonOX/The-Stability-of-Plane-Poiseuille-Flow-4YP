% Computers the operators to be used in the inequality

function [A1,A2,E,W,I] = compute_operators(M,N) %,y,D2,D4,I,W,gamma,tau)

% Create differentiation matricies
[D,y,D2,D4,I,W] = compute_cheb_diff_matricies(N);

% Large matrix to store the matricies along its diagonal
A1 = zeros((N-1)*(2*M+1), (N-1)*(2*M+1));
A2 = zeros((N-1)*(2*M+1), (N-1)*(2*M+1));
E = zeros((N-1)*(2*M+1), (N-1)*(2*M+1));

for n = -M:1:M
    
    % Counter used to store values
    j = n + M + 1;
        
    A1(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = ...
        (-n^4*I + 2*n^2*D2 - D4);   
    
    A2(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = ...
        (2i*I*n + 1i*diag(1 - y(2:N).^2)*(D2*n - I*n^3)); 
    
    E(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = ...
        (I*n^2 - D2);
            
end

% Create WW matrix 
W = kron(eye(2*M+1),W);
I = eye(size(W,1));

% A1 = sparse(A1);
% A2 = sparse(A2);
% E = sparse(E);
% W = sparse(W);
% I = sparse(I);

end