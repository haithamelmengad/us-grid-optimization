function [cost, ghg, var] = calcImpacts(x)
% FUNCTION CALCIMPACTS
% For a given x vector containing values for the 61 decision
% variables, return the respective expected cost,
% GHG emissions (MT), and variance.

load all_problem_data.mat

civVar = civSD .^ 2; 
cicVar = cicSD .^ 2;
ckdVar = ckdSD .^ 2;

%% Set up Cost Calculation

f1 = NaN(61,1);

% First 54 for x_it
ctr = 1; % set counter to index f vector
for i = 1:I
    for t = 1:T
        f1(ctr) = civBar(i) * n_t(t);
        ctr = ctr + 1;
    end
end

% Next 4 for y_i
for i = 5:8
    f1(ctr) = 1000 * cicBar(i);
    ctr = ctr + 1;
end

% last 3 for z_k
for k = 1:K
    f1(ctr) = ckdBar(k) * (sMax_kt(k,:) * n_t);
    ctr = ctr + 1;
end

%% Set up GHG Vector
f2 = NaN(61,1);

% First 54 for x_it
ctr = 1; % set counter to index f vector
for i = 1:I
    for t = 1:T
        f2(ctr) = 0.001 * g_i(i) * n_t(t);
        ctr = ctr + 1;
    end
end

% Next 4 for y_i
for i = 5:8
    f2(ctr) = 0;
    ctr = ctr + 1;
end

% last 3 for z_k
for k = 1:K
    f2(ctr) = 0;
    ctr = ctr + 1;
end

%% Set up Variance Matrix
% From 4.33
% set H according to the fomulation
H3 = zeros(16,61);

% First 9 rows for xit

ctr = 1; 
for i = 1:I
    for t = 1:T
        H3(i,ctr) = n_t(t);
        ctr = ctr + 1;
    end
end

% Next 4 for y_i

for i = 5:8
    H3(9+i-4,ctr) = 1000;
    ctr = ctr + 1;
end

% last 3 for z_k

for k = 1:K
    H3(9+4+k,ctr) = sMax_kt(k,:) * n_t;
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
Hvar = 2 * H3' * sigma * H3;

%% Do the calculations
cost = f1' * x;
ghg = f2' * x;
var = .5 * x' * Hvar * x;

