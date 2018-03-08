% SC, HEM, KX
% 498 SIS Project Part 4.31
% Add Constraints of Emission Cap

% First load the problem setup data
load all_problem_data

%% ----------------------------------------------------%
% 1. How should the cap be set in order to achieve the same solution
% as the minimum environmental emissions objective?
% -----------------------------------------------------%

% The cap should be set as the fval of pt429: MT662825 CO2
% when the cap is set smaller than it, there will no solution
% because the constraint requires the function to be lower than
% its minimun value

% when the cap is set larger than it, the solution will not change

%% ----------------------------------------------------%
% 2. under what values would the emissions cap have an influence
% in the minimum expected cost optimization problem?
% -----------------------------------------------------%

% get the f from 4.29
f = NaN(61,1);

% First 54 for x_it
ctr = 1; % set counter to index f vector
for i = 1:I
    for t = 1:T
        f(ctr) = 0.001 * g_i(i) * n_t(t);
        ctr = ctr + 1;
    end
end

% Next 4 for y_i
for i = 5:8
    f(ctr) = 0;
    ctr = ctr + 1;
end

% last 3 for z_k
for k = 1:K
    f(ctr) = 0;
    ctr = ctr + 1;
end

% set emission cap constraints as 0.001 * gi * nt * xit <= b_cap 
% the constraints are same as the objective function of 429

% set the new constraint by adding a line of f' to A
A_plusCap = [A;f'];


% set new f according to pt428 and get the x and fval of min cost problem
% get the result from pt428
f = NaN(61,1);

% First 54 for x_it
ctr = 1; % set counter to index f vector
for i = 1:I
    for t = 1:T
        f(ctr) = civBar(i) * n_t(t);
        ctr = ctr + 1;
    end
end

% Next 4 for y_i
for i = 5:8
    f(ctr) = 1000 * cicBar(i);
    ctr = ctr + 1;
end

% last 3 for z_k
for k = 1:K
    f(ctr) = ckdBar(k) * (sMax_kt(k,:) * n_t);
    ctr = ctr + 1;
end

%--------------------------------------------------------------
%-----run linear code to get the results of 428---------
%--------------------------------------------------------------
% Set up options: use simplex
options = optimoptions('linprog', 'Algorithm', 'dual-simplex');

% Run different versions depending on Matlab version
verStr = version;
verS = strtok(verStr);
ver = str2double(verS(1:3));
    
if ver >= 9.1 % version 2016
    [x, fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, lb, ub, options);
else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, lb, ub, [], options);
    % Basically adds a x0 vector input before options
end


% ---------------------------------------------------------------------
% ---run linear code with the new constraint on min emission problem---
% ---------------------------------------------------------------------
N = 10^3; % N represents the simulation interval
n = 601; % n represents the numbers of intervals
count = 1; 

% Initialize storage matrix for simulation values
% Col 1: cap constraint, changes each iteration
% Col 2: fval result from linprog
% Col 3: fval under cap - fval min cost
f_cap = NaN(n,3); % [cap, new fval, difference]

for i = 500000:N:1100000 % try different gap of cap for simulation
    cap = i; % define the cap
    b_plusCap = [b; cap]; % add cap constraint to the new b

    if ver >= 9 % version 2016
        [x_plusCap, fval_plusCap, exitflag_plusCap, output_plusCap, lambda_plusCap] = ...
            linprog(f, A_plusCap, b_plusCap, Aeq, beq, lb, ub, options);
    else % for Ketong's 2014. assuming anything before 2016 is like this
        [x_plusCap, fval_plusCap, exitflag_plusCap, output_plusCap, lambda_plusCap] = ...
            linprog(f, A_plusCap, b_plusCap, Aeq, beq, lb, ub, [], options);
        % Basically adds a x0 vector input before options
    end

    % record results
    if exitflag_plusCap == 1
    f_cap(count,:) = [cap, fval_plusCap, fval_plusCap-fval];
    else
        f_cap(count,:) = [cap, NaN, NaN];
    end
    count = count + 1;
end

%% plot the results, cap and the fval
x1=650000:N:1100000;
y1=5.5031e+07 + 0 * x1;

figure(8)
pl = plot(f_cap(:,1),f_cap(:,2),'.',x1,y1,'--');
set(pl,{'markers'}, {20; 20});
xlabel('Emission Cap (MT CO2e)', 'FontSize', 14)
ylabel('Expected Cost ($)', 'FontSize', 14)
title('Emission Cap Embedded in Cost Optimization Problem', ...
    'FontSize', 16, 'fontWeight', 'bold')
hl = legend('Expected Cost', 'Minimum Cost');
set(hl, 'FontSize', 14);
% the result shows that when the cap is set under 1028000
% the minimun expected cost optimization problem would be influnced
