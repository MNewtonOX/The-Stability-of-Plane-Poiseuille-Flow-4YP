% Compares the Reynolds number across different discretization grids
clear all

%% Define Parameters

% Set N values
N_start = 2;
N_step = 1;
N_end = 20;

% Set M values
M_start = 1;
M_step = 1;
M_end = 20;

% Define matrix of Reynolds nummbers
Re = zeros((round((N_end - N_start)/N_step) + 1),(round((M_end - M_start)/M_step) + 1));

for N = N_start:N_step:N_end
    
for M = M_start:M_step:M_end

countN = (N - N_start)/N_step + 1;
countM = (M - M_start)/M_step + 1;
    
%% Call CVX package in semi-positive definite mode

cvx_begin sdp

%% Define varaibles
 
variable lamda(1,1) nonnegative
variable eta(1,1) nonnegative

%% Create Q matrix, such that psi' * Q * psi < 0

[A1,A2,W] = compute_kron_matricies(N,M);

%% Minimise lamda and set constraints

minimise( lamda )

subject to 

    % Constraint in this form to prevent errors in the way CVX works
    W*(lamda*A1 + A2) + (lamda*A1 + A2)'*W' + eta*eye(size(W)) <= 0
    
cvx_end

% Calculate the reynolds number as reciprocal of lamda
Re(countN,countM) = 1/cvx_optval;

end

end
 
