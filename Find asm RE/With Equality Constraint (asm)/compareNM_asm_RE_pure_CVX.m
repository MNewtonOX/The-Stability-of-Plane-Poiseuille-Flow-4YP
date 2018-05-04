% Finds the Reynolds number where the flow is asymptotically stable via the 
% popov cirerion. Computes over varrying discretization grid.
% The answer should be the same as the circle criterion Reynolds nummber.
clear all
cvx_solver SDPT3
cvx_save_prefs

%% Define Parameters
% Set accuracy (number of iterations)
accuracy = 6;

for N = 2:7
    
for M = 1:7

% Assume Re is between 70 and 130
lower_Re = 70;
upper_Re = 120;

for m = 1:accuracy
     
% Try Re for the upper bound
Re = 0.5*(lower_Re + upper_Re);
gamma = 1/Re;

% Create differentiation matricies
[A1,A2,E,W,I] = compute_operators(M,N);

% Compute matrix G to find Gk
[G] = compute_equality_constraint_matrix(M,N);

%% Define CVX variables
cvx_begin sdp

variable Q((2*M+1)*(N-1),(2*M+1)*(N-1)) hermitian
variable tau(1,1) nonnegative
variable eta1(1,1) nonnegative % Eta used to reduce errors with inequality constraint
variable eta2(1,1) nonnegative

minimise ( eta1 + eta2 ) 

subject to

%% First condition (not needed)
L1((M*(N-1) + 1):(3*M+1)*(N-1), (M*(N-1) + 1):(3*M+1)*(N-1)) = ...
(E'*W*Q*E) + (E'*W*Q*E)';
[L1, zeros(size(L1,1),1); zeros(1,size(L1,1)), 0] - eta1*ones(size(L1,1)+1) >= 0; 

%% Second condition
%L = zeros(2*(4*M+1)*(N-1),2*(4*M+1)*(N-1));
WW = kron(eye(2),(W*(-2*I)));
WW = WW(1:(4*M+1)*(N-1), 1:(4*M+1)*(N-1));

% Top left
L2((M*(N-1) + 1):(3*M+1)*(N-1), (M*(N-1) + 1):(3*M+1)*(N-1)) = ...
    (gamma*A1'*W*Q*E + A2'*W*Q*E + gamma*E'*W*Q*A1 + E'*W*Q*A2);

% Top right
L2((M*(N-1) + 1):(3*M+1)*(N-1), ((5*M+1)*(N-1) + 1):(7*M+2)*(N-1)) = ...
    (E'*W*Q - gamma*A1'*W - A2'*W - tau*W*I);

% Bottom left
L2(((5*M+1)*(N-1) + 1):(7*M+2)*(N-1), (M*(N-1) + 1):(3*M+1)*(N-1)) = ...
    (W*Q*E - gamma*W*A1 - W*A2 - tau*W*I);

% Bottom right
L2(((4*M+1)*(N-1) + 1):(2*(4*M+1)*(N-1)), ((4*M+1)*(N-1) + 1):(2*(4*M+1)*(N-1))) = ...
    WW;

[L2, zeros(size(L2,1),1); zeros(1,size(L2,1)), 0] + ...
[L2, zeros(size(L2,1),1); zeros(1,size(L2,1)), 0]' + ...
eta2*ones(size(L2,1) + 1) <= 0;


%% Third condition
for k = 1:(4*M+1)*(N-1)
    
    % Rearrange equations to get into quadratic form
    % Since psi'*Gk*psi + ek*u1 = 0 then x'*Pk*x + hk*x = 0
    
    % Encase Gk in zeros
    Gk = [zeros((M*(N-1)),((4*M+1)*(N-1))); ...
          zeros(((2*M+1)*(N-1)),((M)*(N-1))), G(:,:,k), zeros(((2*M+1)*(N-1)),((M)*(N-1)));...
          zeros((M*(N-1)),((4*M+1)*(N-1)))];
    
    Pk = [Gk, zeros(((4*M+1)*(N-1)),((4*M+1)*(N-1))); ...
         zeros(((4*M+1)*(N-1)),((4*M+1)*(N-1))), zeros(((4*M+1)*(N-1)),((4*M+1)*(N-1)))];
    ek = [zeros(k-1,1); 1; zeros((4*M+1)*(N-1)-k,1)];
    hk = [zeros((4*M+1)*(N-1),1);ek];
    
    % Use S-procedure to create matrix equality
    [Pk 0.5*hk; 0.5*hk' 0] == 0;
        
end

cvx_end

if isnan(cvx_optval) == 1 || cvx_optval == Inf  %cvx_status == 'Solved'
    
    upper_Re = 0.5*(lower_Re + upper_Re);
    
else
    
    lower_Re = 0.5*(lower_Re + upper_Re);

end
clear L1 L2

end

ReS(N,M) = lower_Re;

end

end

