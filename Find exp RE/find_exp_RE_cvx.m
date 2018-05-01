% Calculates the maximum Reynolds number for which flow is expontialy
% stable using a Lyapunov energy method. 
% It uses CVX along with either SDPT3, mosek or SeDuMi. 
clear all

%% Define Parameters
% Set discritisation grid
N = 10; % y direction
M = 10; % x direction

%% Call CVX package in semi-positive definite mode
cvx_begin sdp

%% Define varaibles
 
variable lamda(1,1) nonnegative
variable eta(1,1) nonnegative

%% Create L matrix, such that psi' * L * psi < 0

[A1,A2,W] = compute_kron_matricies(N,M);

%% Minimise lam%da and set constraints
minimise( lamda )

subject to 

    % Constraint in this form to prevent errors in the way CVX works
    W*(lamda*A1 + A2) + (lamda*A1 + A2)'*W' + eta*ones(size(W)) <= 0
    
cvx_end

% Calculate the reynolds number as reciprocal of lamda
Re = 1/cvx_optval;

disp(['Value of Reynolds number is ', num2str(Re)]);
 
