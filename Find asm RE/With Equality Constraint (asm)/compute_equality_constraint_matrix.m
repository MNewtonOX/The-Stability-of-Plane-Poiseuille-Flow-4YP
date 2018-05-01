%% Create the equality constraint matrix
% A ones/zeros (J) matrix is created and used to select which derivative
% terms are to be used from each operator
% Heavily brute force method and is not what you would call elegant

function [G] = compute_equality_constraint_matrix(M,N);

%% Create the ones and zeros matrix
% Preallocate J Matrix
J = zeros((2*M+1)*(N-1),(2*M+1)*(N-1),(4*M+1)*(N-1));

% Create vector of x,y values corresponding to the u vector to be
% referenced later
k = 1;
for i = -2*M:1:2*M
    for j = 1:N-1
        
        u_xy(k,1) = i;
        u_xy(k,2) = j;
        k = k + 1;
        
    end
end

% Go through all elements in u1
for k = 1:(4*M+1)*(N-1)
    
    ux_term = u_xy(k,1);    % Get x reference
    uy_term = u_xy(k,2);    % Get y reference

    m = 1;  % Dummy count variable
    for i = -M:1:M    % Cycle through all combinations of g and h vectors
        for j = -M:1:M
        
            if i + j == ux_term     % Find combinations that sum to k
                
                gx_terms(m) = i + M + 1;    % Save the combination of terms
                hx_terms(m) = j + M + 1;    % as a reference
                m = m + 1;
                
            end
            
        end        
    end
    
    % Populate the J matrix
    for i = 1:size(gx_terms,2)
        
        row_pos = (1 + (gx_terms(i) - 1)*(N-1)) + uy_term - 1; % x part + y part
        col_pos = (1 + (hx_terms(i) - 1)*(N-1)) + uy_term - 1;
        
        J(row_pos,col_pos,k) = 1;
        
    end
    clear gx_terms; % Clear term combinations to avoid issue with next iteration
    clear hx_terms;
end

%% Create differentiation matrix

% Compute standard diff matricies
[D,y,D2,D4,I,W] = compute_cheb_diff_matricies(N);
D3 = D^3;
D = D(2:N,2:N);
D3 = D3(2:N,2:N);

for n = -M:1:M
    
    % Counter used to store values
    j = n + M + 1;
        
    % Creates the operator matricies which are reference to as B    
    B1(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = (1i*n*I);   
    B2(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = (D);
    B3(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = (D3 - n^2*D);
    B4(((j-1)*(N-1) + 1):(j*(N-1)), ((j-1)*(N-1) + 1):(j*(N-1))) = (-1i*n^3*I + 1i*n*D2);
                
end

%% Use ones and zeros and diff matrix to populate G matrix

% Pre-allocate G matrix
G = zeros((2*M+1)*(N-1),(2*M+1)*(N-1),(2*M+1)*(N-1));

for k = 1:(4*M+1)*(N-1)

    % Multiple them together to get linear combination of terms, conjugate
    % of first term must be used to incorporate the right terms
    G(:,:,k) = B1.'*J(:,:,k)*B3 - B2.'*J(:,:,k)*B4; 
    
end


