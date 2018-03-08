% SC, HEM, KX
% 498 SIS Project Part 4.33
% Quadratic Programming for min variance of the annual cost

% first load the problem setup data
load all_problem_data

% work in terms of variances
civVar = civSD .^ 2; 
cicVar = cicSD .^ 2;
ckdVar = ckdSD .^ 2;

%% totCostVar = genCostVar + buildCostVar + DSMCostVar;

% when civBar is multipied by MWh, should multiply
% civVar by the the square of MWh
% genCostVar = (nt.^2) * (xit.^2) * civVar;
% buildCostVar = cicVar(5:8)' * (yi.^2) * (1000^2);
% DSMCostVar = ckdVar' * (zk.^ 2) * (sMax_kt(k,:) .^2 * nt.^2);
% 2 * (H'*sigma*H] + c' * H
% Matlab quadprog form: 0.5 * x' * H * x + f' * x
% there are only quadratic term in the objective function so f' = 0
f = zeros(1,61)';

% set H according to the fomulation
H = zeros(16,61);

% First 9 rows for xit

ctr = 1; 
for i = 1:I
    for t = 1:T
        H(i,ctr) = n_t(t);
        ctr = ctr + 1;
    end
end

% Next 4 for y_i

for i = 5:8
    H(9+i-4,ctr) = 1000;
    ctr = ctr + 1;
end

% last 3 for z_k

for k = 1:K
    H(9+4+k,ctr) = sMax_kt(k,:) * n_t;
    ctr = ctr + 1;
end

sigma = zeros(16);

ctr = 1; % set counter to index f vector
for i = 1:9
    sigma(ctr,ctr) = civVar(i);
    ctr = ctr + 1;
end

for i = 5:8
    sigma(ctr,ctr) = cicVar(i);
    ctr = ctr + 1;
end

for k = 1:3
    sigma(ctr,ctr) = ckdVar(k);
    ctr = ctr +1;
end

% This is the full matrix going into the quadprog
% which is HT E H
H1 = 2 * H' * sigma * H;


%% --------------------
% ---------------------
% Quadratic Programming CODE
% ---------------------
% ---------------------

% Use quadprog to solve min variance function
% Set up options
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
[x, fval, exitflag, output, lambda] = ...
    quadprog(H1, f, A, b, Aeq, beq, lb, ub, [], options);

% result: the min variance is 2.0035e+13
