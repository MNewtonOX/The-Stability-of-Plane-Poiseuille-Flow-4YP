% Tests the Reynolds number where the flow is exponetially stable via the 
% popov cirerion, without using the equality constraint 
% The answer should be the same as the circle criterion Reynolds nummber.
clear all

%% Define Parameters
% Set discritisation grid
M = 2; % x direction
N = 10; % y direction

% Set initial Reynolds number and gamma
Re = 87.5; %85.3
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


% % First condition
% (W*(E'*Q*E)) + (W*(E'*Q*E))' - eta1*ones(size(W)) >= 0;
% 
% % Second condition
% [ (W*(gamma*A1'*Q*E + A2'*Q*E + gamma*E'*Q*A1 + E'*Q*A2))...
%   (W*(E'*Q - gamma*A1' - A2' - tau*I));...
%   (W*(Q*E - gamma*A1 - A2 - tau*I))  ...         
%   (W*(-2*I)) ] + ...
% [ (W*(gamma*A1'*Q*E + A2'*Q*E + gamma*E'*Q*A1 + E'*Q*A2))...
%   (W*(E'*Q - gamma*A1' - A2' - tau*I));...
%   (W*(Q*E - gamma*A1 - A2 - tau*I))  ...         
%   (W*(-2*I)) ]'+ eta2*ones(2*size(W)) <= 0;

