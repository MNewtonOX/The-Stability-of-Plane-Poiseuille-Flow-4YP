% Script to calculate the value of Reynolds number (Re) that causes the
% system to become unstable

% Specify order of accuracy of solution
order = 6;

% Dummy variables
calc_Re = zeros(order,1);
current_Re = 0;

for m = 1:order
    
    for Re = (current_Re + 10^(4-m)):10^(4-m):(current_Re + 10^(5-m))

        % Finds the maximum eigenvalue for a given Re
        [max_eig] = find_max_eigenvalue_n1(Re);

        if max_eig >= 0 
            
            calc_Re(m) = Re - 10^(4-m);
            current_Re = calc_Re(m);
            break

        end
        
    end

end

Re = current_Re;
disp(Re);


