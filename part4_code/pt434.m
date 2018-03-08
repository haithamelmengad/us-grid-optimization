% SC, HEM, KX
% 498 SIS Project Part 4.34
% Multi-objective (Cost variance and Expected cost) Optimization Problem

% first load the problem setup data
load all_problem_data

% work in terms of variances
civVar = civSD .^ 2; 
cicVar = cicSD .^ 2;
ckdVar = ckdSD .^ 2;

%% Find:
% min (theta * f1(x) + (1 - theta) * f2(x))
% Where:
% f1 and f2 are the cost and variance functions

% Theta is the weighting factor between f1 and f2
% Vary theta from 0 to 1 to create Pareto frontier
% theta = 0:N:1

% see Pt 4.28 & 4.33 for optimization problem setup

% get H1 and f1 for cost, see pt428
H1 = zeros(61);
f = NaN(1,61);

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

f1 = f';

% get H2, f2, see pt433

f2 = zeros(61,1);

H = zeros(16,61);
% First 9 rows for xit

ctr = 1; % set counter to index f vector
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

zigma = zeros(16);

ctr = 1; 
for i = 1:9
    zigma(ctr,ctr) = civVar(i);
    ctr = ctr + 1;
end

for i = 5:8
    zigma(ctr,ctr) = cicVar(i);
    ctr = ctr + 1;
end

for k = 1:3
    zigma(ctr,ctr) = ckdVar(k);
    ctr = ctr +1;
end

H2 = 2 * H' * zigma * H;

% N represent the step size of theta
N = 0.01;
numSteps = length(0:N:1);

% Initialize a cell array to store quadprog results
multObjRes_434 = cell(6, numSteps);

% Storage matrix for f1' and f2' results
% column 1: theta
% column 2: 1-theta
% column 3: fval - cost
% column 4: fval - variance
fComp = NaN(numSteps,4);
count = 1;

%% --------------------------------------------------
% for loop to change theta to obtain different x and fval
% then record fval1 and fval2
% --------------------------------------------------


for theta = 0:N:1
    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
    [x, fval, exitflag, output, lambda] = ...
        quadprog(theta*10^-7*H1+(1-theta)*10^-13*H2, theta*10^-7*f1+(1-theta)*10^-13*f2, A, b, Aeq, beq, lb, ub, [], options);

    % Save values
    multObjRes_434{count, 1} = x;
    multObjRes_434{count, 2} = fval;
    multObjRes_434{count, 3} = exitflag;
    multObjRes_434{count, 4} = output;
    multObjRes_434{count, 5} = lambda;
    multObjRes_434{count, 6} = theta;
    
    % Store both fvals in a matrix for charting
    fComp(count,:) = [theta, 1-theta, 0.5*x'*H1*x+f1'*x, 0.5*x'*H2*x+f2'*x];
    count = count + 1;
end

%% ---------------------
% ---------------------
% PLOT PARETO FRONTIER
% ---------------------
% ---------------------
figure(1)
plot(fComp(:,3),fComp(:,4),'-k');
xlabel('f1(x): Cost')
ylabel('f2(x): Variance')
title('Pareto Frontier of Efficient Solutions to Min Cost and Variance (Normalized)')

hold on
scatter(fComp(:,3),fComp(:,4),20,fComp(:,1),'filled'); 
colorbar
ylabel(colorbar,'Theta')

% the graph shows that the fluctuate points are clustering in
% one consecutive line. The theta goes from 0.32 to 0.47;
% the cost goes from 62743980.21 to 62361907.54
% the variance goes from 2.34996E+13 to 2.37878E+13

% record the x from the medium results (theta = 0.39) as x_pt434
x_pt434 = multObjRes_434{40,1};