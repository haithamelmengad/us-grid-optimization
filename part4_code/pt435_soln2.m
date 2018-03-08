% Dec 12 - after fixing variance
% Looks like the tradeoffs between 2 of the 3
% in the triple objective is much more apparent now.


% load data from pt435
load 435_fvalComp
load 435_tripleObjRes

% Optimal results from 4.28, 4.29, and 4.33 respectively
minCost = 5.5031E07;
minGHG = 6.6283E05;
minVar = 2.0035E13;

%% Screen fvalComp
% Check if there are any values in fvalComp
% that are within x% of all 3 min objectives.

tol = 1.35; % tolerance multiplier

% Make 3 bool vectors to screen fvalComp for each criteria
boolCost = fvalComp(:, 4) <= tol*minCost;
boolGHG = fvalComp(:, 5) <= tol*minGHG; 
boolVar = fvalComp(:, 6) <= tol*minVar;

% 
% See how many meet all 3
meetAll = boolCost + boolGHG + boolVar;
boolAll = meetAll == 3;

% History:
% First try at 25%. OK totals at 25% over were 3058, 2249, 827.
% but sum boolAll is zero. max meetAll was 2. Need to expand.
% Try 2: 30%. sum boolAll is still zero.
% Try 3: 35%. sum boolAll now 18.

% Save these solutions.
solSet = fvalComp(boolAll == 1, :);
solSetOptResults = tripleObjRes(boolAll == 1, :);

% Plot these 18 solutions at 35% tolerance
figure(20)
scatter3(solSet(:,4), solSet(:,5), solSet(:,6))
rotate3d on 
xlabel('Expected Cost')
ylabel('GHG')
zlabel('Cost Variance')
% OK they are pretty scattered, no clear trend.
% Do other analyses from pt435_soln.m - first run

%% Find min reliance on renewables (total MW)

% Renewables are i = 7, 8, 9
% Corresponding x vector indices are 37-54

numSolsTest = size(solSet, 1);

% storage matrices
renewableMW = NaN(numSolsTest, 1);

for i = 1:numSolsTest
    % extract the x vector
    testxvec = solSetOptResults{i, 1};
    renSum = sum(testxvec(37:54));
    renewableMW(i) = renSum;
end

% Do a little histogram
figure(21)
histogram(renewableMW, 20)
% fairly spread in small range

% min = 1.5909E03
% max = 1.6720E03
% mean = 1.6269E03




%% Find Min New Constructions
% Calculate the total new construction
% indices 55-58 are greater than zero, then count it
% Also sum it for total size of new constructions


% storage matrices
buildCount = NaN(numSolsTest, 1);
newMWSum = NaN(numSolsTest, 1);

for i = 1:numSolsTest
    % Extract the y values
    testxvec = solSetOptResults{i, 1};
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
figure(22)
histogram(newMWSum,20)
% spread, but somewhat bottom heavy.

% min = 265.1
% max = 357.7
% mean = 304
% less than 285: 7


%% Selection and Analysis
% Do a chart comparing newgen and renewables min
figure(24)
scatter(newMWSum, renewableMW, 'filled')
xlabel('Total New Capacity Constructed (MW)')
ylabel('Amount of Renewable Energy Used (MW)')
title('Selected Solutions: Renewables Used vs. New Capacity Constructed')

% OK there is a pretty linear line, to be expected
% since most of the renewables are new construction

% select the bottom left point: lowest new capacity
% and 3rd lowest renewables used
% Index: 18 (in this 18-value matrix)
% New Capacity: 265 MW
% Renewables Used: 1.593 E03 MW

% Extract corresponding x vector and calculate impacts
indSel = 18; % selected index
xSel = solSetOptResults{indSel, 1};
[selCost, selGHG, selVar] = calcImpacts(xSel);
% selCost = 6.6414E07
% selGHG = 8.8492E05
% selVar = 2.6976E13

% Find percent differences
diffCost = (selCost-minCost)/minCost;
diffGHG = (selGHG-minGHG)/minGHG;
diffVar = (selVar-minVar)/minVar;
% 20.7%, 33.5%, 34.7%



