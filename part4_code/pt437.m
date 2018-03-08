% pt437.m
% Checking sensitivity by increasing ou_i(9) from 8% to 62%
% CEE 498 SIS Project
% S Cai, K Xie, H El Mengad

load all_problem_data
load 435_x_newdesign
[cost_design, ghg_design, var_design] = calcImpacts(x_newdesign);

%% Change constraints related to the ou(9) and find the new design
% first get the new ou_i vector
new_ou_i = ou_i;
new_ou_i(9) = 0.62;

% the change of ou_i(9) only affect the upper bounds of instant capacity
new_ub = ub;
new_ub(46:54) = (1 - new_ou_i(9)) * xMax_i(9);

% run through part 435 again to obtain the new solution

%----------------------------------------------------------%
%--------------------------pt435---------------------------%
H1 = zeros(61);
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

H2 = zeros(61);

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

% work in terms of variances
civVar = civSD .^ 2; 
cicVar = cicSD .^ 2;
ckdVar = ckdSD .^ 2;

% No linear decision variables, so f vector is zeros
f3 = zeros(61,1);

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

Hvar = 2 * H3' * sigma * H3;
c1 = 1E-7;
c2 = 1E-5;
c3 = 1E-13;

% Define the three scaling factors: theta, phi, and 1-theta-phi
% Store values in a 3-column matrix
stepSize = 0.01;
numRuns = 5151; % based on an initial run to find number
scaleFactors = NaN(numRuns, 3); % Store factor values
ctr = 1; % counter to go through scaleFactors storage matrix

for theta = 0:stepSize:1
    for phi = 0:stepSize:1
        % Check to make sure theta + phi is not > 1
        if (theta+phi) > 1
            continue
        end
        % find 3rd factor
        remainder = 1-theta-phi;
        scaleFactors(ctr, 1) = theta;
        scaleFactors(ctr, 2) = phi;
        scaleFactors(ctr, 3) = remainder;
        ctr = ctr + 1;
    end
end

tripleObjResHydro = cell(numRuns, 5);

fvalCompHydro = NaN(numRuns,5);
for i = 1:numRuns
    % Grab scaling factor values
    fac1 = scaleFactors(i, 1);
    fac2 = scaleFactors(i, 2);
    fac3 = scaleFactors(i, 3);
    
    % Set up equations
    HTemp = fac1.*c1.*H1 + fac2.*c2.*H2 + fac3.*c3.*Hvar;
    fTemp = fac1.*c1.*f1 + fac2.*c2.*f2 + fac3.*c3.*f3;
    
    % Run optimization
    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
    [x, fval, exitflag, output, lambda] = ...
        quadprog(HTemp, fTemp, A, b, Aeq, beq, lb, new_ub, [], options);

    % Save values
    tripleObjResHydro{i, 1} = x;
    tripleObjResHydro{i, 2} = fval;
    tripleObjResHydro{i, 3} = exitflag;
    tripleObjResHydro{i, 4} = output;
    tripleObjResHydro{i, 5} = lambda;
    
    % Store results in charting matrix
    fvalCompHydro(i,1) = fac1;
    fvalCompHydro(i,2) = fac2;
    fvalCompHydro(i,3) = fac3;
    % These are individual calcs
    fvalCompHydro(i,4) = (0.5*x'*H1*x) + (f1'*x);
    fvalCompHydro(i,5) = (0.5*x'*H2*x) + (f2'*x);
    fvalCompHydro(i,6) = (0.5*x'*Hvar*x) + (f3'*x);
end
%--------------------------pt435---------------------------%
%----------------------------------------------------------%

% run through pt435_soln2.m
%------------------------------------------------------------%
%----------------------pt435.soln2.m-------------------------%
% results:
% tol should be set to 1.5
% There are six points

% Optimal results from 4.28, 4.29, and 4.33 respectively
minCost = 5.5031E07;
minGHG = 6.6283E05;
minVar = 2.0035E13;

tolHydro = 1.5; % tolerance multiplier, there are six solutions

% Make 3 bool vectors to screen fvalComp for each criteria
boolCostHydro = fvalCompHydro(:, 4) <= tolHydro*minCost;
boolGHGHydro = fvalCompHydro(:, 5) <= tolHydro*minGHG; 
boolVarHydro = fvalCompHydro(:, 6) <= tolHydro*minVar;

% 
% See how many meet all 3
meetAllHydro = boolCostHydro + boolGHGHydro + boolVarHydro;
boolAllHydro = meetAllHydro == 3;

% Save these solutions.
solSetHydro = fvalCompHydro(boolAllHydro == 1, :);
solSetOptResHydro = tripleObjResHydro(boolAllHydro == 1, :);

% Find min reliance on renewables (total MW)

% Renewables are i = 7, 8, 9
% Corresponding x vector indices are 37-54

numSolsTestHydro = size(solSetHydro, 1);

% storage matrices
renewableMWHydro = NaN(numSolsTestHydro, 1);

for i = 1:numSolsTestHydro
    % extract the x vector
    testxvec = solSetOptResHydro{i, 1};
    renSum = sum(testxvec(37:54));
    renewableMWHydro(i) = renSum;
end


% Find Min New Constructions

% storage matrices
buildCountHydro = NaN(numSolsTestHydro, 1);
newMWSumHydro = NaN(numSolsTestHydro, 1);

for i = 1:numSolsTestHydro
    % Extract the y values
    testxvec = solSetOptResHydro{i, 1};
    yVals = testxvec(55:58);
    % Find number of new builds
    numBuilds = sum(yVals > 0);
    % Find total new MW
    totNewMW = sum(yVals);
    % Store values
    buildCountHydro(i) = numBuilds;
    newMWSumHydro(i) = totNewMW;
end

% Selection and Analysis
% Do a chart comparing newgen and renewables min
figure(24)
scatter(newMWSumHydro, renewableMWHydro, 'filled')
xlabel('Total New Capacity Constructed (MW)')
ylabel('Amount of Renewable Energy Used (MW)')
title('Selected Solutions: Renewables Used vs. New Capacity Constructed')

% Index: 2 or 4

% Extract corresponding x vector and calculate impacts (x-61 and
% cost,emissions and variance)

impactHydro = NaN(61+3, numSolsTestHydro);
for indSel = 1:numSolsTestHydro % selected index
    xSelHydro = solSetOptResHydro{indSel, 1};
    [selCostHydro, selGHGHydro, selVarHydro] = calcImpacts(xSelHydro);
    impactHydro(:,indSel) = [xSelHydro;selCostHydro;selGHGHydro;selVarHydro];
end

% see excel for the results


%% Add constraints of the new construction based on our new design
% Based on our construction decision, i.e. y_i results
% Use equality constraints to force y_i values
yi = x_newdesign(55:58);
new_Aeq = zeros(4,61);
new_Aeq(1,55) = 1;
new_Aeq(2,56) = 1;
new_Aeq(3,57) = 1;
new_Aeq(4,58) = 1;
new_beq = yi;

% run through part 435 again to obtain the new solution

tripleObjResHydroCon = cell(numRuns, 5);

fvalCompHydroCon = NaN(numRuns,5);
for i = 1:numRuns
    % Grab scaling factor values
    fac1 = scaleFactors(i, 1);
    fac2 = scaleFactors(i, 2);
    fac3 = scaleFactors(i, 3);
    
    % Set up equations
    HTemp = fac1.*c1.*H1 + fac2.*c2.*H2 + fac3.*c3.*Hvar;
    fTemp = fac1.*c1.*f1 + fac2.*c2.*f2 + fac3.*c3.*f3;
    
    % Run optimization
    [x, fval, exitflag, output, lambda] = ...
        quadprog(HTemp, fTemp, A, b, new_Aeq, new_beq, lb, new_ub, [], options);

    % Save values
    tripleObjResHydroCon{i, 1} = x;
    tripleObjResHydroCon{i, 2} = fval;
    tripleObjResHydroCon{i, 3} = exitflag;
    tripleObjResHydroCon{i, 4} = output;
    tripleObjResHydroCon{i, 5} = lambda;
    
    % Store results in charting matrix
    fvalCompHydroCon(i,1) = fac1;
    fvalCompHydroCon(i,2) = fac2;
    fvalCompHydroCon(i,3) = fac3;
    % These are individual calcs
    fvalCompHydroCon(i,4) = (0.5*x'*H1*x) + (f1'*x);
    fvalCompHydroCon(i,5) = (0.5*x'*H2*x) + (f2'*x);
    fvalCompHydroCon(i,6) = (0.5*x'*Hvar*x) + (f3'*x);
end

% run through pt435_soln2.m
% tol shoule be set to 1.52

tolHydroCon = 1.52; % tolerance multiplier

% Make 3 bool vectors to screen fvalComp for each criteria
boolCostHydroCon = fvalCompHydroCon(:, 4) <= tolHydroCon*minCost;
boolGHGHydroCon = fvalCompHydroCon(:, 5) <= tolHydroCon*minGHG; 
boolVarCon = fvalCompHydroCon(:, 6) <= tolHydroCon*minVar;

% 
% See how many meet all 3
meetAllHydroCon = boolCostHydroCon + boolGHGHydroCon + boolVarCon;
boolAllHydroCon = meetAllHydroCon == 3;

% Save these solutions.
solSetCon = fvalCompHydroCon(boolAllHydroCon == 1, :);
solSetOptResultsCon = tripleObjResHydroCon(boolAllHydroCon == 1, :);
numSolsTestHydroCon = size(solSetCon, 1);

% Extract corresponding x vector and calculate impacts (x-61 and
% cost,emissions and variance)

impactHydroCon = NaN(61+3, numSolsTestHydroCon);
for indSel = 1:numSolsTestHydroCon % selected index
    xSelHydroCon = solSetOptResultsCon{indSel, 1};
    [selCostHydroCon, selGHGHydroCon, selVarHydroCon] = ...
        calcImpacts(xSelHydroCon);
    impactHydroCon(:,indSel) = ...
        [xSelHydroCon;selCostHydroCon;selGHGHydroCon;selVarHydroCon];
end

% Combining all the results together
impactTotal = [[x_newdesign;cost_design;ghg_design;var_design],...
    impactHydro,impactHydroCon];
% see excel for the results
