% SC, HEM, KX
% 498 SIS Project Part 4.28
% Minimum Expected Cost Optimization Problem

% Find:
% a Min cost
% b Which plants & DSM programs used?
% c CO2eq emissions associated with design?
% d SD of cost?
% e Which constraints satisfied at equality?
% f Sensitivity analysis - 
%   find constraint with highest shadow price

% See formulation.m for optimization problem setup and load
% all_problem_data.mat to get the data
load all_problem_data

%% --------------------
% ---------------------
% OBJECTIVE FUNCTION f
% ---------------------
% ---------------------

% Decision variables
% x_it = power produced from plant i in load block t (MW)
% y_i = design capacity of new plant i
% z_k = implementation rate of DSM program k
% 54, 4, and 3 DV's respectively = 61 total
% Min linear function, use linprog
% Get into form: f' x
% Get all the decision variables into a column vector
% x = [x11, x12, ... x16, x21, x22 ... x26 ... x96,
%   y5, y6, y7, y8, z1, z2, z3]
% Index for y5: 55
% Index for z1: 59

% Now set up f transpose
% first 54 values of f' are cvi * nt: cv1 * n1, cv1 * n2,... 
% cv1*n6, cv2*n1, ... cv2*n6... up to cv9 * n6
% next 4 are 1000 * cc[5,6,7,8]
% last 3 are c1*(n1s11 + n2s12 + ... + n6s16) etc

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

%% --------------------
% ---------------------
% LINEAR OPTIMIZATION CODE
% ---------------------
% ---------------------

% Use linprog to solve min cost function
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

% fval result: the min cost is $5.5031e+07

% Group data into an array and export
% SYC 12/05
minExpCostResults = cell(1, 5);
minExpCostResults{1,1} = x;
minExpCostResults{1,2} = fval;
minExpCostResults{1,3} = exitflag;
minExpCostResults{1,4} = output;
minExpCostResults{1,5} = lambda;
save('minExpCostResults.mat', 'minExpCostResults');


%% --------------------
% ---------------------
% CO2 EMISSIONS
% ---------------------
% ---------------------

% trim x solution to first 54 values to get x_it only
xitOpt = x(1:54);
% make it a 9x6 matrix
xitOptM = reshape(xitOpt, [6,9])';
% multiply by hours per load block to get MWh
xitOptMwh = xitOptM * n_t;
% multiply by g_i * 0.001 for tons of CO2e
minCostGHG = xitOptMwh' * g_i * 0.001;
% result: 1.0279e+06


%% --------------------
% ---------------------
% COST STANDARD DEVIATION
% ---------------------
% ---------------------
% work in terms of variances
civVar = civSD .^ 2; 
cicVar = cicSD .^ 2;
ckdVar = ckdSD .^ 2;

% When values are scaled, variances are scaled by the square
% then you can sum them together for total variance
% Assuming costs are SIID
% ie when civBar is multipied by MWh, should multiply
% civVar by the the square of MWh
genCostVar = civVar' * (xitOptMwh.^2);
genCostSD = sqrt(genCostVar);
% result: 6.1816e+06

% y_i
% extract values from x
yiOpt = x(55:58);
% strip down variance vector to 4 values
cicVar = cicVar(5:8);
% multiply out
buildCostVar = cicVar' * (yiOpt.^2) * 1000^2;
buildCostSD = sqrt(buildCostVar);
% result: 8.0968e+05

% z_k
% extract values from x
zkOpt = x(59:end);
% multiply out
DSMCostVar = ckdVar' * (zkOpt .^ 2) * (sMax_kt(k,:) * n_t)^2;
DSMCostSD = sqrt(DSMCostVar); 
% result: 1.7000e+5

% Add variances together for total cost SD
totCostVar = genCostVar + buildCostVar + DSMCostVar;
totCostSD = sqrt(totCostVar);
% result: 6.2367e+06
% Basically latter two are inconsequential in total

%% --------------------
% ---------------------
% Find Slack Variables
% ---------------------
% ---------------------

% Find which constraints have slack variables
% by testing the optimal solution by all the constraints (39+122)
% and calculating the difference
% There are 39 constraints in A matrix
% and 61x2 through the upper & lower bounds
slack = NaN(39+122, 1);
ctr = 1;
% find the slack for each constraint (39+122)
for i = 1:39
    slack(ctr) = b(ctr) - A(ctr,:) * x;
    ctr = ctr + 1;
end

for j = 1:61
    slack(ctr) = x(ctr - 39) - lb(ctr - 39);
    ctr = ctr + 1;
end

for k = 1:61
    slack(ctr) = ub(ctr - 39 - 61) - x(ctr - 39 - 61);
    ctr = ctr + 1;
end

% find the constrainsts that are satisfied at equally
bind_c = abs(slack)<=10^-14;
sum(bind_c)

% find the constraints having slack variables
s_v = abs(slack)>10^-14;
sum(s_v)

% results: there are 76 constrainsts that are satisfied at equally

%% --------------------
% ---------------------
% Sensitivity Analysis Code
% Calculate Shadow Price
% ---------------------
% ---------------------

% For each constraint (39+122 in total), calculate shadow price:
% change of fval / change of the constraint

shadow_price = NaN(39+122, 1);
n1 = 0.01; % n1 = change of th constraint, because zk varies from
           % 0~1, therefore the change of the lower bound constraints
           % should be set smaller than 1 and as small as possible.
           % This doen't affect other constraints
ctr = 1;

% Use for loop to calculate each shadow price for each constraint (39+122)
for i = 1:39 % the first 39 constraints
    new_b = b;
    new_b(i) = new_b(i) + n1; % new b under the change of ith constraint
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex');
    % Version of Matlab for linprog
    if ver >= 9.1 % version 2016
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, new_b, Aeq, beq, lb, ub, options);
    else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, new_b, Aeq, beq, lb, ub, [], options);
    % Basically adds a x0 vector input before options
    end
    n2 = new_fval - fval; % n2 = change of fval
    shadow_price(ctr) = n2/n1;
    ctr = ctr + 1;
end

for j = 1:61 % following 61 lower bound constraints
    new_lb = lb;
    new_lb(j) = new_lb(j) + n1; % new lb from the change of jth constraint
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex');
    % Version of Matlab for linprog
    if ver >= 9.1 % version 2016
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, new_lb, ub, options);
    else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, new_lb, ub, [], options);
    % Basically adds a x0 vector input before options
    end
    n2 = new_fval - fval; % n2 = change of fval
    shadow_price(ctr) = n2/n1;
    ctr = ctr + 1;
end

for k = 1:61 % following 61 upper bound constraints
    new_ub = ub;
    new_ub(k) = new_ub(k) + n1; % new ub from the change of kth constraint
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex');
    % Version of Matlab for linprog
    if ver >= 9.1 % version 2016
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, lb, new_ub, options);
    else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, lb, new_ub, [], options);
    % Basically adds a x0 vector input before options
    end
    n2 = new_fval - fval; % n2 = change of fval
    shadow_price(ctr) = n2/n1;
    ctr = ctr + 1;
end
% Find which constraint has the highest shadow price
find(shadow_price==max(shadow_price))

% Result: the 99th constraint has the highest sensitivity;
% it is the 60th lb constraint: z_2 >= 0

%% Combine the results
sensitivity = NaN(161,2);
sensitivity(:,1) = bind_c;
sensitivity(:,2) = shadow_price;
