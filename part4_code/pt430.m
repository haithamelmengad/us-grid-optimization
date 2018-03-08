% SC, HEM, KX
% 498 SIS Project Part 4.30
% Multi-objective (Emissions and Expected cost) Optimization Problem

% Find:
% min (theta * f1(x) + (1 - theta) * f2(x))
% Where:
% f1 and f2 are the emission and cost functions
% % From 4.21 f1 = emission
% % f1 = 0.001 * gi * nt * xit
% % From 4.22 f2 = Cost
% % f2 = (cvi * nt * xit + 1000 * cc[5,6,7,8] * yi) + ...
% % ck * (n1s11 + n2s12 + ... + n6s16) * zk

% Theta is the weighting factor between f1 and f2
% Vary theta from 0 to 1 to create Pareto frontier
% theta = 0:N:1

% see Pt 4.28 & 4.29 for optimization problem setup
load all_problem_data

%% get f1, see pt429
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

f1 = f;

% get f2, see pt428
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

f2 = f;

% N represent the step size of theta
N = 0.001;
numSteps = length(0:N:1);

% Initialize a cell array to store linprog results
multObjRes_430 = cell(6, numSteps);

% Storage matrix for f1' and f2' results
fComp = NaN(numSteps,4);
count = 1;

%% -------------------------------------------------
% for loop to change theta to obtain different x and fval
% then record each f1*x and f2*x
% --------------------------------------------------
options = optimoptions('linprog', 'Algorithm', 'dual-simplex');

for theta = 0:N:1


% Run different versions depending on Matlab version
verStr = version;
verS = strtok(verStr);
ver = str2double(verS(1:3));

if ver >= 9.1 % version 2016
    [x, fval, exitflag, output, lambda] = ...
        linprog(theta*f1*10^-5+(1-theta)*f2*10^-7, A, b, Aeq, beq, lb, ub, options);
else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, fval, exitflag, output, lambda] = ...
        linprog(theta*f1*10^-5+(1-theta)*f2*10^-7, A, b, Aeq, beq, lb, ub, [], options);
    % Basically adds a x0 vector input before options
end

    % Save values
    multObjRes_430{count, 1} = x;
    multObjRes_430{count, 2} = fval;
    multObjRes_430{count, 3} = exitflag;
    multObjRes_430{count, 4} = output;
    multObjRes_430{count, 5} = lambda;
    multObjRes_430{count, 6} = theta;
    
    % Store both fvals in a matrix for charting
    fComp(count,:) = [theta, 1-theta, f1' * x, f2' * x];
    count = count + 1;
end

%% --------------------
% ---------------------
% PLOT PARETO FRONTIER
% ---------------------
% ---------------------
figure(1)
plot(fComp(:,3),fComp(:,4),'-k');
xlabel('Global Warming Potential (MT CO2e)', ...
    'FontSize', 14)
ylabel('Expected Cost ($)', 'FontSize', 14)
title('Pareto Frontier of Efficient Solutions (Normalized)'...
    , 'FontSize', 16, 'fontWeight', 'bold')
hold on
scatter(fComp(:,3),fComp(:,4),40,fComp(:,1),'filled'); 
colorbar
ylabel(colorbar,'Theta')

% The result shows a clear fluctuate, and save the x corresponding to
% the medium one as x_pt430 
x_pt430 = multObjRes_430{823,1};