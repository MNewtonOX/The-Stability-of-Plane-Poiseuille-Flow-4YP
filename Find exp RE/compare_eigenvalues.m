% Show eigenvalues are all negative across different grid sizes
clear all;

% True Reynolds number
Re = 87.7;
lamda = 1/Re;

N_test = 10:5:40;

j = 1;
% Testing over different grid sizes
for N = N_test
    
    M = N ;

    % Compute matricies
    [A1,A2,W] = compute_kron_matricies(N,M);
    
    Q = W*(lamda*A1 + A2) + (lamda*A1 + A2)'*W';
    
    % Only take the maximium eigenvalue as that will determine stability
    e(j) = max(eig(Q));
    
    j = j + 1;
    
end

% Plot graph of eigenvalues
plot(N_test',e);

