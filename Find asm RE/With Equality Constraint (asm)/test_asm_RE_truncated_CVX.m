% Tests the Reynolds number where the flow is asymptotically stable via the 
% popov cirerion, but by truncating the u1 vector. Computes over varrying
% discretization grid.
% The answer should be the same as the circle criterion Reynolds nummber.
clear all

%% Define Parameters
% Set discritisation grid
M = 3; % x direction
N = 8; % y direction

% Set initial Reynolds number and gamma
Re = 86.25; %85.3
gamma = 1/Re;

% Create differentiation matricies
[A1,A2,E,W,I] = compute_operators(M,N);

% Compute matrix G to find Gk
[G] = compute_equality_constraint_matrix(M,N);

% Truncate G to get the relevant terms
G = G(:,:,(M*(N-1) + 1):((3*M+1)*(N-1)));

%% Define CVX variables
cvx_begin sdp

variable Q((2*M+1)*(N-1),(2*M+1)*(N-1)) hermitian
variable tau(1,1) nonnegative
variable eta1(1,1) nonnegative % Eta used to reduce errors with inequality constraint
variable eta2(1,1) nonnegative

minimise ( eta1 + eta2 ) 

subject to

%% First condition
L1 = (E'*W*Q*E) + (E'*W*Q*E)';
[L1, zeros(size(L1,1),1); zeros(1,size(L1,1)), 0] - eta1*ones(size(L1,1) + 1) >= 0;

%% Second condition
L2 = [ (gamma*A1'*W*Q*E + A2'*W*Q*E + gamma*E'*W*Q*A1 + E'*W*Q*A2),...
  (E'*W*Q - gamma*A1'*W - A2'*W - tau*W*I);...
  (W*Q*E - gamma*W*A1 - W*A2 - tau*W*I), ...         
  (-2*W*I) ] + ...
[ (gamma*A1'*W*Q*E + A2'*W*Q*E + gamma*E'*W*Q*A1 + E'*W*Q*A2),...
  (E'*W*Q - gamma*A1'*W - A2'*W - tau*W*I);...
  (W*Q*E - gamma*W*A1 - W*A2 - tau*W*I), ...         
  (-2*W*I) ]';

[L2, zeros(size(L2,1),1); zeros(1,size(L2,1)), 0] + eta2*ones(size(L2,1) + 1) <= 0;

%% Third condition
for k = 1:(2*M+1)*(N-1)
    
    % Rearrange equations to get into quadratic form
    % Since psi'*Gk*psi + ek*u1 = 0 then x'*Pk*x + hk*x = 0
    Pk = [G(:,:,k) zeros(size(G(:,:,k))); zeros(size(G(:,:,k))) zeros(size(G(:,:,k)))];
    ek = [zeros(k-1,1); 1; zeros((2*M+1)*(N-1)-k,1)];
    hk = [zeros((2*M+1)*(N-1),1);ek];
    
    % Use S-procedure to create matrix equality
    S = sparse([Pk 0.5*hk; 0.5*hk' 0]);
    S == 0;
        
end

cvx_end
