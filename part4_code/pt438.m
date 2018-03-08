% SC, HEM, KX
% 498 SIS Project Part 4.30
% Multi-objective (Emissions and Expected cost) Optimization Problem

% Find:
% min (theta * f1(x) + (1 - theta) * f2(x))
% Where:
% f1 and f2 are the emission and cost functions
% % From 4.21 f1 = emission
% % f1 = 0.001 * gi * nt * xit
% % From 4.22 f2 = Cost
% % f2 = (cvi * nt * xit + 1000 * cc[5,6,7,8] * yi) + ...
% % ck * (n1s11 + n2s12 + ... + n6s16) * zk
% theta = 0.4

%---------------------------------------------------------------------
%              SETUP
%---------------------------------------------------------------------

% load all_problem_data for problem set up
load all_problem_data

% according to the prompt, there are two values changed in the original 
% formulation due to the techonology break through of CCS
g_i(5) = 16;
cicBar(5) = 43.4;

%% get f1, see pt429
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

f1 = f;

% get f2, see pt428
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

f2 = f;

%% -------------------------------------------------
% Multi-objective optimization on cost and emissions
% using linprog
% --------------------------------------------------
options = optimoptions('linprog', 'Algorithm', 'dual-simplex');
count = 1; %initialize counter
N=0.001;

for theta = 0:N:1


% Run different versions depending on Matlab version
verStr = version;
verS = strtok(verStr);
ver = str2double(verS(1:3));

if ver >= 9.1 % version 2016
    [x, fval, exitflag, output, lambda] = ...
        linprog(theta*f1*10^-5+(1-theta)*f2*10^-7, A, b, Aeq, beq, lb, ub, options);
else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, fval, exitflag, output, lambda] = ...
        linprog(theta*f1*10^-5+(1-theta)*f2*10^-7, A, b, Aeq, beq, lb, ub, [], options);
    % Basically adds a x0 vector input before options
end

    % Save values
    multObjRes_438{count, 1} = x;
    multObjRes_438{count, 2} = fval;
    multObjRes_438{count, 3} = exitflag;
    multObjRes_438{count, 4} = output;
    multObjRes_438{count, 5} = lambda;
    multObjRes_438{count, 6} = theta;
    
    count = count + 1;
end

% Store results for optimized cost and emissions for theta = 0.6 --
% 60% emissions 40% cost weighting scheme
numVector = [multObjRes_438{:,6}]; %creates a vector based on the 6th column of the multObjRes_438 array
index = find(numVector==0.4); %indexes the element of numVector equal to 0.6 (theta)
x_438 = multObjRes_438{(index),1}; %assigns the x value that corresponds to theta=0.6 to x_438
result_438 = [ f1' * x_438, f2' * x_438]; %total cost and total emissions associated with the result

%---------------------------------------------------------------------
%                ANALYSIS
%---------------------------------------------------------------------

%Analyze Advanced Coal with CCS (ACCS) usage:

x_438_comp = [ multObjRes_438{:,1}]; %store x vectors for different thetas in a matrix called x_438_comp
CCSnew_use = sum(x_438_comp(25:30,:)); %Stores the total ACCS use in each scenario in vector CCSnew_use
%column 1 is total ACCCS use for theta= 0.20, column 7 is ACCS use for theta= 0.80,
CCSnew_pos = find(CCSnew_use>0); %stores the elements that are superior to 0
size = numel(CCSnew_pos); %sets size equal to the number of elements in CCSnew_pos
thetas_CCS = NaN(size,1); %creates an empty matrix the size of CCnew_pos
count= 1 ; %initiate counter
for i = 1001-size:1000
thetas_CCS(count,:) = [multObjRes_438{i,6}]; %gives corresponding theta values for which CCS use is superior to 0.
count = count + 1;
end
plot(thetas_CCS*100, CCSnew_pos) %plot Advanced Coal with CCS use versus theta values
xlabel('Percent weighting in favor of Emissions') %x axis title
ylabel ('Total Advanced Coal with CCS usage in MW') %y axis title
title ('Impact on the change of theta on Advanced Coal with CCS usage, after CCS technology improvement') %chart title

%------------------------------------------------------------------
%               COMPARISON WITH pt4.30 RESULTS
%------------------------------------------------------------------

%% extract the results from problem 4.30 with theta = 0.60
load 430_multObjRes.mat
numVector = [multObjRes_430{:,6}]; %creates a vector based on the 6th column of the multObjRes_430 array
index = find(numVector==0.6); %indexes the element equal to 0.6 (theta)
x_430 = multObjRes_430{(index),1}; %assigns the x value that corresponds to theta to x_430
%results are presented in Table 438, Advanced Coal with CCS use is 0
%% extract the results from problem 4.30 with theta = 0, 0.10, 0.20,0.50, 0.80, 0.90 and 1
numVector = [multObjRes_430{:,6}]; %creates a vector based on the 6th column of the multObjRes_430 array
index = [find(numVector==0), find(numVector==0.1), find(numVector==0.2), find(numVector==0.5), find(numVector==0.8), find(numVector==0.9), find(numVector==1) ]; %indexes the elements equal to 0.2, 0.5, 0.7(thetas)
x_430_comp = [ multObjRes_430{index(1),1} multObjRes_430{index(2),1} multObjRes_430{index(3),1} multObjRes_430{index(4),1} multObjRes_430{index(5),1} multObjRes_430{index(6),1} multObjRes_430{index(7),1}]; %stores the associated x solutions in Matrix form
%column 1 is x solution for theta= 0.20, column 2 is x solution for theta= 0.50,
%column 3 is x solution for theta = 0.70
CCS_use = sum(x_430_comp(25:30,:)); %gives the total Adv. Coal CCS use for each of the three scenarios
%interpretation: Adv. Coal isn't used before tech. advancements whether emissions are prioritized or not


