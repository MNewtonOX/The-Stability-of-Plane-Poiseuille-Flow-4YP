% Tests if the equality constraint matrix is valid by comparing to know
% solutions for the stream functions and u1.
clear all

% Set grid
M = 3;
N = 5;

% Create x and y points
x = (-M*2*pi):(2*pi):(M*2*pi); % Are not used, just for reference
[~,y] = cheb(N); 

%% Calculate analytically
% The calculations for this are in my log book and are needed to understand
% the problem

% Create psi vector for all x,y
% psi = cos(3x)*(1-y^2)^2
psi = zeros((2*M+1)*(N-1),1);
u1a = zeros((4*M+1)*(N-1),1);

for i = -2*M:1:2*M
    for j = 2:(length(y) - 1) % Removes boundary conditions
    
        if i == 6
        
            % Calculate row position
            row = (i+2*M)*(N-1) + (j - 1);

            % u1a calculated by hand using derivatives
            u1a(row,1) = -9i*y(j)*(5*(y(j)^2) + 3)*(1 - (y(j)^2))^3;     
                       
        elseif i == -6
                                   
            row = (i+2*M)*(N-1) + (j - 1);

            u1a(row,1) = 9i*y(j)*(5*(y(j)^2) + 3)*(1 - (y(j)^2))^3;  
            
        elseif i == -3 
            
            row = (i+M)*(N-1) + (j - 1);
            
            % psi that fits the BC
            psi(row,1) = 0.5*(1 - y(j)^2)^3;  
            
        elseif i == 3
            
            row = (i+M)*(N-1) + (j - 1);
            
            % psi that fits the BC
            psi(row,1) = 0.5*(1 - y(j)^2)^3;   
            
        else
            
            row = (i+2*M)*(N-1) + (j - 1);

            u1a(row,1) = 0; 

        end

    end
end


%% Calculate with equality constraint matrix

% Get equality matrix
[G] = compute_equality_constraint_matrix(M,N);

% Truncate G to get the relevant terms
%G = G(:,:,(M*(N-1) + 1):((3*M+1)*(N-1)));

for k = 1:(4*M+1)*(N-1)

    % Compute using psi from before
    u1b(k,1) = -psi'*G(:,:,k)*psi;
        
end

%% Compare methods

% Difference between two vectors
diff = u1a - u1b;

% Correlation of vectors (normalized mse)
R = corrcoef(u1a(abs(u1a) > 0.00001),u1b(abs(u1b) > 0.00001));

% Fit the u1b to a polynomial to compare coefficents
p = polyfit(y(2:N),u1b(1:N-1),9);
p(2,1:9) = abs(-9:-1);


