% Script to calculate the value of Reynolds number (Re) that causes the
% system to become unstable by looking at all the nodes

% Specify order of accuracy of solution
order = 4;

% Dummy variables
calc_Re = zeros(order,1);
current_Re = 0;

for m = 1:order
    
    for Re = (current_Re + 10^(4-m)):10^(4-m):(current_Re + 10^(5-m))

        % Finds the maximum eigenvalue for a given Re
        [max_eig,n] = find_max_eigenvalue_n_vary_0_2(Re);

        if max_eig > 0   %unsure if can also equal zero
            
            calc_Re(m) = Re - 10^(4-m);
            current_Re = calc_Re(m);
            break

        end
        
    end

end

% Display and store the value of Re
Re = current_Re;
text1 = ['Value of Reynolds number is ', num2str(Re)];
text2 = ['Becomes unstable when n = ', num2str(n)];
disp(text1);
disp(text2);

