% pt435_soln.m
% Following the triple-objective optimzation in pt435.m
% Calculate the corresponding selected solution.

% assume matrix 'fvalComp' and cell 'tripleObjRes'
% from 4.35 optimization is in memory

% 2. Select final values and find a x-vector
% that corresponds to final values.

% Target cost ~ $8E07 and 7E5 GHG
targetCost = 8E7;
targetGHG = 7E5;

numRuns = 5151; % empirical from last part

%% Extract set of close solutions
% First extract values from corresponding to these two
% and all associated variances
% Use about 10% values in either direction

margin = 0.1;


% OLD code that just found a list of values
% What I actually need are indices
% so I can extract the x values from tripleObjRes
% targetCostVals = NaN(1, 6);
% ctr = 1;
% for i = 1:size(fvalComp,1)
%     testCost = fvalComp(i, 4);
%     if (testCost >= (1-margin)*targetCost) && ...
%         (testCost <= (1+margin)*targetCost)
%             targetCostVals(ctr,:) = fvalComp(i, :);
%             ctr = ctr+1;
%     end
% end
% 
% 
% % Now strip out everything that doesn't meet
% % GHG bounds
% delRows = zeros(size(targetCostVals, 1),1);
% for i = 1:size(targetCostVals,1)
%     testGHG = targetCostVals(i, 5);
%     if (testGHG <= (1-margin)*targetGHG) || ...
%             (testGHG >= (1+margin)*targetGHG)
%         delRows(i) = 1; % delete
%     end
% end
% 
% targetCGVals = targetCostVals(delRows == 0, :);
% 
% % 283 remaining variances to compare


% New version:
% Find indices where costs are within range
costLo = fvalComp(:, 4) >= (1-margin)*targetCost;
costHi = fvalComp(:, 4) <= (1+margin)*targetCost;

% Find indices where GHG in range
ghgLo = fvalComp(:, 5) >= (1-margin)*targetGHG;
ghgHi = fvalComp(:, 5) <= (1+margin)*targetGHG;

% Combine
allGood = costLo + costHi + ghgLo + ghgHi;
% values are 4 if everything matches
% extract into binary vector for simplicity
allGoodBin = allGood == 4;

% size is 283 after just controlling for cost and GHG
% Now find cost and variance of each, and plot


%% Select Optimal Variance Values

% Subset of fvalComp with cost & GHG within 10%
allGoodVals = fvalComp(allGoodBin == 1, :);

% Do a quick plot
figure(8)
histogram(allGoodVals(:, 6), 20)

% OK, weird distribution, gonna get rid of <1E13

% Subset with cost & GHG within 10%, and lower vars
allGoodVarVals = allGoodVals(allGoodVals(:, 6) <= ...
    1E13, :);
% Now left with 111

% plot this histogram
figure(9)
histogram(allGoodVarVals(:, 6), 20)

% Find the corresponding binary vector
% With length 5151 to mask the full results
allGoodVarBin = fvalComp(:, 6) <= 1E13;
allGood2 = allGood + allGoodVarBin;
allGoodBin2 = allGood2 == 5;
% sum is 111, all good

% OK let's go with the min new constructions
% as a heuristic, then select the min variance of those
% first I need to extract all the x_its


%% Find Min New Constructions

% Extract the corresponding 111 x vectors
% Using the binary mask vectory allGoodBin2
scrnSols = tripleObjRes(allGoodBin2 == 1, :);

% -------------------------------
% Calculate the total new construction
% indices 55-58 are greater than zero, then count it
% Also sum it for total size of new constructions
numSolsTest = size(scrnSols, 1);

% storage matrices
buildCount = NaN(numSolsTest, 1);
newMWSum = NaN(numSolsTest, 1);

for i = 1:numSolsTest
    % Extract the y values
    testxvec = scrnSols{i, 1};
    yVals = testxvec(55:58);
    % Find number of new builds
    numBuilds = sum(yVals > 0);
    % Find total new MW
    totNewMW = sum(yVals);
    % Store values
    buildCount(i) = numBuilds;
    newMWSum(i) = totNewMW;
end

% OK, everything requires new constructions
% IE all of buildCount is == 4

% Do a quick histogram of new MW
figure(10)
histogram(newMWSum,20)
% sum(newMWSum <= 450) = 3
% sum(newMWSum <= 440) s= 1

%% Find min reliance on renewables (total MW)
% Renewables are i = 7, 8, 9
% Corresponding x vector indices are 37-54

% storage matrices
renewableMW = NaN(numSolsTest, 1);

for i = 1:numSolsTest
    % extract the x vector
    testxvec = scrnSols{i, 1};
    renSum = sum(testxvec(37:54));
    renewableMW(i) = renSum;
end

% Do a little histogram
figure(11)
histogram(renewableMW, 20)

% sum(renewableMW <= 1550) = 21
% sum(renewableMW <= 1540) = 9
% sum(renewableMW <= 1530) = 5

% Do a chart comparing newgen and renewables min
figure(12)
scatter(newMWSum, renewableMW)
xlabel('Total New Capacity Constructed (MW)')
ylabel('Amount of Renewable Energy Used (MW)')

% OK there is a pretty straight line now... 

[minNewMW, mnmInd] = min(newMWSum);
% 114.9, 1

[minRenMW, mrmInd] = min(renewableMW);
% 1.0360E03, 1

selSolInd = 1;

%% Selected Solution Calculations

% Extract the optimal solution
xSel = scrnSols{selSolInd, 1};
fvalSel = scrnSols{selSolInd, 2}; % 2.8662
flagSel = scrnSols{selSolInd, 3}; % 1

% Calculate expected cost, GHG, and variance
% Assume f1, f2, and H3 from 4.35 are in memory
selCost = f1' * xSel; % 7.5908 E07 $, 37% increase
selGHG = f2' * xSel; % 1.1720 E06 MT, 77% increase
selVar = 0.5 * xSel' * H3 * xSel; % 2.8662 E12 $^2





