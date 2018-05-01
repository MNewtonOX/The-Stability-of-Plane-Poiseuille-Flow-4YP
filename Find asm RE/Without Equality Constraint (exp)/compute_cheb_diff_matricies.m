% Creates the differentiation matricies to save space and repeating the
% same code

function [D,y,D2,D4,I,W] = compute_cheb_diff_matricies(N)

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
W = diag(w(2:N));

end