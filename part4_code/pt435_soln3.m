%% Adjusting the pt435 solution to be more realistic

% from the result we discover that there are some unrealistic
% situation like building 0.28 MW capcity solar plant
% therefore, we decide to change it to zero while checking if
% the new solution still meet all the constraints by calculating
% the slack variables

% first load the data from previous solution and problem setup
load 435_xSel
load all_problem_data

minCost = 5.5031E07;
minGHG = 6.6283E05;
minVar = 2.0035E13;

selCost = 6.6414E07;
selGHG = 8.8492E05;
selVar = 2.6976E13;

% adjusting the xit and yi that less than 1 to be zeros
% (basically from the solar plant)
x_adjust = xSel(1:58);
x_adjust(x_adjust <= 1) = 0;
x_adjust = [x_adjust; xSel(59:61)];

% using the code from pt428/429 to check the slack variables
% in order to find out if the new solution meet all the constraints
slack = NaN(39+122, 1);
ctr = 1;
% find the slack for each constraint (39+122)
for i = 1:39
    slack(ctr) = b(ctr) - A(ctr,:) * x_adjust;
    ctr = ctr + 1;
end

for j = 1:61
    slack(ctr) = x_adjust(ctr - 39) - lb(ctr - 39);
    ctr = ctr + 1;
end

for k = 1:61
    slack(ctr) = ub(ctr - 39 - 61) - x_adjust(ctr - 39 - 61);
    ctr = ctr + 1;
end

% the result shows that the new solution doesn't meet the first four
% load blocks due to we omitted solar energy which is basically used
% in these four load blocks.
% the slack variables are quite small (0.28 at maximum)
% so we decide to use advance CC2 natural gas to fill the gaps
% while not exceeding the upper bound

% fill first-third load block
x_adjust(31:33) = x_adjust(31:33) + 0.28;
% fill fourth load block
x_adjust(34) = x_adjust(34) + 0.17;
% Change the new y6 capacity
x_adjust(56) = 93.96/(1-0.05);

% Then test again the slack variables and they are all positive now
slack = NaN(39+122, 1);
ctr = 1;
for i = 1:39
    slack(ctr) = b(ctr) - A(ctr,:) * x_adjust;
    ctr = ctr + 1;
end

for j = 1:61
    slack(ctr) = x_adjust(ctr - 39) - lb(ctr - 39);
    ctr = ctr + 1;
end

for k = 1:61
    slack(ctr) = ub(ctr - 39 - 61) - x_adjust(ctr - 39 - 61);
    ctr = ctr + 1;
end

%% Check the cost, emissions and variance of the new solution
[adjcost, adjghg, adjvar] = calcImpacts(x_adjust);

% Find percent differences between the new one and the selected one
diffCost_adj = (selCost-adjcost)/selCost;
diffGHG_adj = (selGHG-adjghg)/selGHG;
diffVar_adj = (selVar-adjvar)/selVar;
% 2.0995e-04, -1.5219e-04 and -9.6112e-05

% Find percent differences
diffCost = (adjcost-minCost)/minCost;
diffGHG = (adjghg-minGHG)/minGHG;
diffVar = (adjvar-minVar)/minVar;
% 20.66%(20.68% before), 33.53%(33.51 before), 34.66%(34.65% before)

% the difference can be accepted, therefore we have our new design
x_newdesign = x_adjust;