% SC, HEM, KX
% 498 SIS Project Part 4.29
% Minimum Environmental Emissions Optimization Problem

% Find:
% a Min Environmental Emissions
% b cost
% c Which plants & DSM programs used?
% d CO2eq emissions associated with design?
% e SD of cost?
% f Which constraints satisfied at equality/have slack?
% g Sensitivity analysis - 
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
% first 54 values of f' are 0.001 * gi * nt: 0.001*g1 * n1, 0.001*g1 * n2,... 
% 0.001*g1*n6, 0.001*g2*n1, ... 0.001*g2*n6... up to 0.001*g9 * n6
% next 4 are 0
% last 3 are 0

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

%% --------------------
% ---------------------
% LINEAR OPTIMIZATION CODE
% ---------------------
% ---------------------

% Use linprog to solve min emission function
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


% fval result: the min emission is MT6.6283e+05 CO2 eq

% Group data into an array and export
% SYC 12/05
minGWPRes = cell(1, 5);
minGWPRes{1,1} = x;
minGWPRes{1,2} = fval;
minGWPRes{1,3} = exitflag;
minGWPRes{1,4} = output;
minGWPRes{1,5} = lambda;
save('minGWPRes.mat', 'minGWPRes');

%% --------------------
% ---------------------
% EXPECTED COST
% ---------------------
% ---------------------

% trim x solution to get x_it, y_i, z_k seperately
xitOpt = x(1:54);
yiOpt = x(55:58);
zkOpt = x(59:61);
% make x_it a 9x6 matrix
xitOptM = reshape(xitOpt, [6,9])';
% multiply by hours per load block to get MWh
xitOptMwh = xitOptM * n_t;
% multiply by cv_i for energy costs
minGHG_GenCost = xitOptMwh' * civBar;
% result: 3.7437e+07

% multiply y_i by cc_i(5:8) for construction costs
minGHG_BuildCost = yiOpt' * 1000 * cicBar(5:8);
% result: 138540000

% multiply z_k by ck_d for DSM costs
minGHG_DSMCost = zkOpt' * ckdBar * (sMax_kt(k,:) * n_t);
% result: 2200000

% add up all the costs
minGHGCost = minGHG_GenCost + minGHG_BuildCost + minGHG_DSMCost;
% result: 1.7818e+08

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
% result: 6.8061e+06

% y_i
% strip down variance vector to 4 values
cicVar = cicVar(5:8);
% multiply out
buildCostVar = cicVar' * (yiOpt.^2) * 1000^2;
buildCostSD = sqrt(buildCostVar);
% result: 1.0306e+07

% z_k
% multiply out
DSMCostVar = ckdVar' * (zkOpt .^ 2) * (sMax_kt(k,:) * n_t)^2;
DSMCostSD = sqrt(DSMCostVar); 
% result: 2.5962e+05

% Add variances together for total cost SD
totCostVar = genCostVar + buildCostVar + DSMCostVar;
totCostSD = sqrt(totCostVar);
% result: 1.2353e+07
% Basically come from the building cost of new plants
 
%% --------------------
% ---------------------
% Find Slack Variables
% ---------------------
% ---------------------
% 
% Find which constraints have slack variables
% by testing the optimal solution by all the constraints (39+122)
% and calculating the difference
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


% results: there are 62 constrainsts that are satisfied at equally
 
 
%% --------------------
% ---------------------
% Sensitivity Analysis Code
% Calculate Shadow Price
% ---------------------
% ---------------------
 
% For each constraint (39+122 in total), calculate shadow price:
% round objective function / round constraint
shadow_price = NaN(39+122, 1);
n1 = 0.01; % n1 = change of constraint, because the zk varies from
          % 0~1, therefore the change of the lower bound constraints
          % should be set smaller than 1 and as small as possible.
          % This doen't affect other constraints
ctr = 1;

% Simulate each shadow price for each constraint (39+122)
for i = 1:39
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
    n2 = new_fval - fval; % n2 = change of objective function
    shadow_price(ctr) = n2/n1;
    ctr = ctr + 1;
end
 
for j = 1:61
    new_lb = lb;
    new_lb(j) = new_lb(j) + n1; % new lb under the change of jth constraint
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
    n2 = new_fval - fval; % n2 = change of objective function
    shadow_price(ctr) = n2/n1;
    ctr = ctr + 1;
end

for k = 1:61
    new_ub = ub;
    new_ub(k) = new_ub(k) + n1; % new lb under the change of jth constraint
    % Version of Matlab for linprog
    if ver >= 9.1 % version 2016
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, lb, new_ub, options);
    else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, new_fval, exitflag, output, lambda] = ...
        linprog(f, A, b, Aeq, beq, lb, new_ub, [], options);
    % Basically adds a x0 vector input before options
    end
    n2 = new_fval - fval; % n2 = change of objective function
    shadow_price(ctr) = n2/n1;
    ctr = ctr + 1;
end

% Find which constraint has the highest sensitivity
find(shadow_price==max(shadow_price))
 
% Result: the 45th constraint has the highest sensitivity;
% it is the 6th lb constraint: x_16 >= xMin_6

%% Combine the results
sensitivity = NaN(161,2);
sensitivity(:,1) = bind_c;
sensitivity(:,2) = shadow_price;