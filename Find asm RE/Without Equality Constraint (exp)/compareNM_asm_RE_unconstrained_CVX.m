% Finds the Reynolds number where the flow is exponetially stable via the 
% popov cirerion, without using the equality constraint with varrying 
% discretization grid.
% The answer should be the same as the circle criterion Reynolds nummber.
clear all
cvx_solver SDPT3
cvx_save_prefs

%% Define Parameters
% Set discritisation grid

% Set accuracy (number of iterations)
accuracy = 8;

for N = 7
    
for M = 1:7

% Assume Re is between 70 and 130
lower_Re = 80;
upper_Re = 120;

for m = 1:accuracy
     
% Try Re for the upper bound
Re = 0.5*(lower_Re + upper_Re);
gamma = 1/Re;

% Create differentiation matricies
[A1,A2,E,W,I] = compute_operators(M,N);

%% CVX
cvx_begin sdp

variable Q((2*M+1)*(N-1),(2*M+1)*(N-1)) hermitian
variable tau(1,1) nonnegative
variable eta1(1,1) nonnegative % Eta used to reduce errors with inequality constraint
variable eta2(1,1) nonnegative

minimise ( eta1 + eta2 ) 

subject to

% First condition
(E'*W*Q*E) + (E'*W*Q*E)' - eta1*ones(size(W)) >= 0;

% Second condition
[ (gamma*A1'*W*Q*E + A2'*W*Q*E + gamma*E'*W*Q*A1 + E'*W*Q*A2),...
  (E'*W*Q - gamma*A1'*W - A2'*W - tau*W*I);...
  (W*Q*E - gamma*W*A1 - W*A2 - tau*W*I), ...         
  (-2*W*I) ] + ...
[ (gamma*A1'*W*Q*E + A2'*W*Q*E + gamma*E'*W*Q*A1 + E'*W*Q*A2),...
  (E'*W*Q - gamma*A1'*W - A2'*W - tau*W*I);...
  (W*Q*E - gamma*W*A1 - W*A2 - tau*W*I), ...         
  (-2*W*I) ]'+ eta2*ones(2*size(W)) <= 0;

cvx_end

if isnan(cvx_optval) == 1 || cvx_optval == Inf  %cvx_status == 'Solved'
    
    upper_Re = 0.5*(lower_Re + upper_Re);
    
else
    
    lower_Re = 0.5*(lower_Re + upper_Re);

end


end

ReS(N,M) = lower_Re;

end

end


