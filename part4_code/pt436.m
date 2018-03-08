% pt436.m
% Calculating final expected cost in $/kW to generate electricity
% CEE 498 SIS Project
% S Cai, K Xie, H El Mengad

load all_problem_data
load 435_x_newdesign
% Calculate the total cost of the new design
[cost, ghg, var] = calcImpacts(x_newdesign);

% Calculate the total capacity of the new design
capacity = xMax_i;
capacity(5:8) = x_newdesign(55:58);
tot_cap = sum(capacity);

% Calculate the total generated electricity
% trim x solution to first 54 values to get x_it only
xit = x_newdesign(1:54);
% make it a 9x6 matrix
xitM = reshape(xit, [6,9])';
% multiply by hours per load block to get MWh
xitMwh = xitM * n_t;

%% Calculate the expected cost
% first way is to calculate by the install capacity
% Result is $/kW
exp_cost1 = cost/tot_cap * 10^-3;

% second way is to calculate by the generated power
% Result is $/kWh
exp_cost2 = cost/sum(xitMwh) * 10^-3;

% The $/kWh excluding capital costs can be calculated as
exp_cost3 = (cost - 1000 * cicBar(5:8)' * x_newdesign(55:58))/...
    sum(xitMwh) * 10^-3;
% the result is 44.12, 0.0193 and 0.0147
