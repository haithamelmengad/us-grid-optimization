% 498 SIS Project Part 2
% SC, HEM, KX

%HEM 11/08 : added documentation
clear all %removes all variables, globals, functions and MEX links
close all %closes all the open figure windows 

% HEM 11/08: deleted the code adding
%the excel file to the working directory
%added documentation

% ======================
% PROBLEM SETUP
% ======================
% Vector x coordinates for each critical infrastructure
% x vector - for reference
% 01 air transp
% 02 electricity
% 03 tlc wireless
% 04 tlc wired
% 05 water mgmt
% 06 rail transp
% 07 finance
% 08 fuel & petroleum grid
% 09 nat gas
% 10 naval ports
% 11 sat comm & nav

% HEM 11/08 : added documentation
% Matrix A - 6-12 Hr Outage Impacts
% Matrix A represents the interdependency coefficients
% A(i,j) is the coefficient of failure propagation from
% infrastructure j to the ith infrastructure
% A=xlsread('Book1.xlsx','Sheet1');
A = [
    0,	0.134,	0.308,	0.456,	0.033,	0.024,	0.012,	0.024,	0.007,	0.001,	0.31;
0,	0,	0.01,	0.023,	0.001,	0.003,	0,	0.002,	0.178,	0.008,	0.004;
0.002,	0.109,	0,	0.12,	0.002,	0.005,	0.002,	0.002,	0.004,	0.002,	0.007;
0.006,	0.083,	0.013,	0,	0.004,	0.003,	0.004,	0.001,	0.004,	0.002,	0.005;
0.005,	0.05,	0.009,	0.02,	0,	0.005,	0.008,	0.008,	0.008,	0.005,	0.02;
0.001,	0.233,	0.104,	0.109,	0.005,	0,	0.007,	0.006,	0.001,	0.002,	0.004;
0.003,	0.1,	0.03,	0.1,	0.007,	0.003,	0,	0.003,	0.007,	0.003,	0.008;
0.008,	0.5,	0.1,	0.05,	0.05,	0.02,	0.02,	0,	0,	0.02,	0.008;
0.002,	0.03,	0.009,	0.005,	0.005,	0,	0.002,	0.005,	0,	0.005,	0.005;
0.005,	0.03,	0.02,	0.02,	0.03,	0.03,	0.02,	0.01,	0.02,	0,	0.005;
0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0;];  %defines Table 2. as Matrix A

save('MatrixA.mat','A'); %saves matrix A file as 'Matrix.mat'



% ======================
% Question 2.8
% ======================
%HEM 11/10 : reviewed documentation changed variable names
% 08 fuel & petroleum grid - find greatest influencer
% find the highest values for row 8
hrow_fuel = find(A(8,:)==max(A(8,:)));
disp(hrow_fuel)
% We observe that sector 2, electricity
% has the greatest influence over the fuel & petroleum grid sector

% 10 naval ports - find greatest influencer
% find the highest values for row 10. 
hrow_ports = find(A(10,:)==max(A(10,:)));
disp(hrow_ports)
% We observe that sectors 2 (elec), 5 (water), and 6 (rail)
% have the greatest influence over naval port sector

% 08 fuel & petroleum grid - find greatest influencee
%find the highest values for column 8. 
hcol_fuel = find(A(:,8)==max(A(:,8)));
disp(hcol_fuel)
% We observe that fuel & petroleum grid
% influences sector 1 (air transportation) the most

% 10 naval ports - find greatest influencee
%find the highest values for column 10.
hcol_ports=find(A(:,10)==max(A(:,10)));
disp(hcol_ports)
% We observe that naval ports influence
% sector 8 (fuel & petroleum grid) the most



% ======================
% Question 2.9
% ======================
% HEM 11/10 : added documentation
% Calculate dependency index kappa and inflence index lambda
% find the dimensions of Matrix A
[m,n] = size(A);
% calculate the two indices of dependency Kappa and influence Lambda
Kappa = sum(A,2)/(m-1); % row summation 
Lambda = sum(A,1)/(n-1); % column summation

% plot Kappa and Lambda
% HEM 11/10 : deleted old line chart code documentation
% Bar chart
% create m-by-2 NaNs matrix for grouped bars
kl = NaN(m,2);
kl(:,1) = Kappa'; %Set 1st column of kl equal to transpose of Kappa
kl(:,2) = Lambda'; %Set 2nd column of kl equal to transpose of Lambda
figure(1)
bar(kl) %create a bar graph of kl
xlabel('Sector')
ylabel('Index Value')
legend('\kappa_{i} Dependency Index',...
    '\lambda_{j} Influence Index')
title('Dependency and Influence Indices for Different Sectors')



% ======================
% Question 2.10
% ======================
% Calculate matrix S = (I-A)^-1 and 
% overall indices kappa_bar, lambda_bar

% create identity matrix
I = eye(n);

% calculate S
S = inv(I-A);
save('MatrixS.mat','S');

% calculate overall indices and plot them
% HEM 11/10 : changed variable names
kappa_overall = sum(S,2)/(m-1); % row summation
lambda_overall = sum(S,1)/(n-1); % column summation 

% HEM 11/10 : deleted old line chart code
% and added documentation

%create a bar plot to compare overall indices to indices
%obtained in Matrix A
klBar = NaN(m,2); %creates a matrix of zeros size m*2 where m is the
                  %number of rows in matrix A
klBar(:,1) = kappa_overall'; %defines kappa_overall transpose as the 1st column of 
                        %klBar
klBar(:,2) = Kappa';%defines Kappa transpose as the 2nd column of
                        %klBar
figure(2)
bar(klBar)  %graphs klBar in a bar plot
xlabel('Sector')
ylabel('Index Value')
leg = legend('$\bar{\kappa}_{i}$ Overall Dependency Index',...
    '$\bar{\lambda}_{j}$ Dependency Index');
set(leg,'Interpreter', 'latex')
title('Overall Dependency Indices vs Dependency Indices for Different Sectors')

klBar = NaN(m,2); %creates a matrix of zeros size m*2 where m is the
                  %number of rows in matrix A
klBar(:,1) = Lambda'; %defines Lambda transpose as the 1st column of 
                        %klBar
klBar(:,2) = lambda_overall';%defines LambdaBar transpose as the 2nd column of
                        %klBar
figure(3)
bar(klBar)  %graphs klBar in a bar plot
xlabel('Sector')
ylabel('Index Value')
leg = legend('$\bar{\kappa}_{i}$ Influence Index',...
    '$\bar{\lambda}_{j}$ Overall Influence Index');
set(leg,'Interpreter', 'latex')
title('Overall Influence Indices vs Influence Indices for Different Sectors')



% ======================
% Question 2.11
% ======================
% HEM 11/10 : reviewed and changed documentation
% Refer to report for justification
% Assume natural gas isn't used for electricity generation to the
% smart grid. Thus dependencies increase for the 3 communications
% infrastructures: 3 TLC wireless, 4 TLC wired, 11 Satellite
% Assume increase of 10% in electricity dependencies on 3,4,11



% ======================
% Question 2.12 - Smart Grid
% ======================
% calculate new A matrix with increased smart grid dependencies
% Assume 10% increase on existing dependencies on 3,4,11

A_sg = A; % smart grid
A_sg(2,3) = 1.1*A(2,3);
A_sg(2,4) = 1.1*A(2,4);
A_sg(2,11) = 1.1*A(2,11);

%calculate new S
S_sg = inv(I-A_sg);

% plot new indices in a bar graph
% HEM 11/10 : changed to bar graph
kappa_new = sum(S_sg,2)/(m-1);
lambda_new = sum(S_sg,1)/(n-1);
kl_new = NaN (m,2);
kl_new(:,1) = kappa_new';
kl_new(:,2) = lambda_new';

figure(4)
bar(kl_new)
xlabel('Sector')
ylabel('Index Value')
leg = legend('$\bar{\kappa}_{i}$ Overall Dependency Index',...
    '$\bar{\lambda}_{j}$ Overall Influence Index');
set(leg,'Interpreter', 'latex')
title('Overall Dependency and Influence Indices for Different Sectors');

% HEM 11/10 : Changed difference plot to bar graph
% Plot the difference between overall indices to interpret
kappa_diff = kappa_new - kappa_overall;
lambda_diff = lambda_new - lambda_overall;
indices_diff = NaN(m,2);
indices_diff(:,1) = kappa_diff;
indices_diff(:,2) = lambda_diff;

figure(5)
bar(indices_diff)
xlabel('Sector')
ylabel('Increase of Overall Indices Value')
legend('Increase of Overall Dependency Index','Increase of Overall Influence Index')
title('Increases of the Two Indices for Different Sectors')



% ======================
% Question 2.13 - Snowstorm
% ======================
% Create vector f_ss to indicate external shock from snowstorm
% 50% damage to 6 rail
% 10% damage to 2 elec
% 35% damage to 1 air
% 20% damage to 10 naval
f_ss = zeros(n,1);
f_ss(6) = 0.5;
f_ss(2) = 0.1;
f_ss(1) = 0.35;
f_ss(10) = 0.2;

% Calculate vector x that indicates infrastructure impacts from snowstorm
x_ss = S * f_ss;



% ======================
% Question 2.15 - Impact Propagation
% ======================
% Simulate propagation of initial damage through kth tier of impacts using
% the recursive relation:
% x(k) = Ax(k-1) + f

% HEM 11/10 : deleted old vector code

k = 10; % number of tiers of impacts to simulate
% Store the x vectors (columns) for each tier
X_ip = NaN(m,k); % NaN for initialization
% X stores corresponding x vector value, NOT correct to store f_ss here
% First column of X_ip is f_ss above
X_ip(:,1) = f_ss;

% Apply recursive function for rest of k
for i = 2:k
    X_ip(:,i) = A * X_ip(:,i-1) + f_ss; 
end

% HEM 11/10 : What's going on here ?
figure(6) %plot x(k) for all the sectors with different color codes
k=0:1:9; % number of kth

subplot(1,2,1)
plot(k,(X_ip(1,:)),'y-')
hold on
plot(k,(X_ip(2,:)),'m-');
hold on
plot(k,(X_ip(6,:)),'y--');
hold on
plot(k,(X_ip(8,:)),'r--');
hold on
plot(k,(X_ip(10,:)),'b--');
xlabel('k')
ylabel('x(k)')
legend('Air Transportation','Elecctricity','Rail Transportation','Fuel&Petro Grid','Naval Ports')

subplot(1,2,2)
plot(k,(X_ip(3,:)),'r-');
hold on
plot(k,(X_ip(4,:)),'g-');
hold on
plot(k,(X_ip(5,:)),'b-');
hold on
plot(k,(X_ip(7,:)),'m--');
hold on
plot(k,(X_ip(9,:)),'g--');
hold on
plot(k,(X_ip(11,:)),'c--');
xlabel('k')
ylabel('x(k)')
legend('TLC Wireless','TLC Wired','Water Management','Finance','Natrual Gas','Satellite Com&Navigation')

% ======================
% Question 2.16 - Uncertainty in Factors
% ======================
% Assume values of A are now random variables
% Uniform distribution +/- 10% of give values
% Use Monte Carlo to determine uncertainty in x given uncertainty in A

% HEM 11/10 : deleted old code

runs = 10000; % number of MC runs
% Initialize storage matrix for x values
xAll = NaN(m,runs);
% Use snowstorm external shock f_ss 
for i = 1:runs
    % Generate matrix of random multipliers
    % Want random values between 0.9 and 1.1, inclusive
    mults = (1.1-0.9) .* rand(m) + 0.9;
    % Use element-wise matrix multiplication to sample a new matrix A
    A_samp = A .* mults; 
    % Calculate x vector 
    x_samp = (I-A_samp)\f_ss;
    % Put into the storage matrix
    xAll(:,i) = x_samp;
end

% HEM 11/10 : added documentation, changed the plot codes to match
% new code
% plot the distributions of uncertainty for each sector in histogram
% with a fit curve to show the probability density function
figure(7)
%suptitle ('Uncertainty distributions of degradation for sectors 1-6')
subplot(3,2,1)
histfit(xAll(1,:));
title('Air Transportation');
subplot(3,2,2)
histfit(xAll(2,:));
title('Electricity');
subplot(3,2,3)
histfit(xAll(3,:));
title('Wireless Telecom');
subplot(3,2,4)
histfit(xAll(4,:));
title('Wired Telecom');
subplot(3,2,5)
histfit(xAll(5,:));
title('Water Management');
subplot(3,2,6)
histfit(xAll(6,:));
title('Rail Transportation');

figure (8)
%suptitle('Uncertainty distributions for sectors 7-11')
subplot (3,2,1)
histfit(xAll(7,:));
title('Finance');
subplot(3,2,2)
histfit(xAll(8,:));
title('Fuel Grid');
subplot(3,2,3)
histfit(xAll(9,:));
title('Natural Gas')
subplot(3,2,4)
histfit(xAll(10,:));
title('Naval Ports')
subplot(3,2,5)
histfit(xAll(11,:));
title('Satellite Com')

% HEM 11/10 : changed variable names
%calculate the mean values
xMean = mean(xAll,2);  
diff = xMean - x_ss;

%find the one sector with largest difference with the expected value
highest_diff = find(abs(diff) == max(abs(diff)));
disp(highest_diff)
% the result shows that the sector rail transportation has the highest diff

%calculate the standard deviation
stdev = std(xAll');
highest_stdev = find((stdev) == max(stdev));
disp(highest_stdev)
% the result shows that the sector fuel&petro grid has the higest std

%Use a box plot to represent uncertainty
figure (9)
subplot(3,2,1)
boxplot(xAll(1,:), 'Air Transportation')
hold on
plot(x_ss(1),'.')
subplot(3,2,2)
boxplot(xAll(2,:), 'Electricity')
hold on
plot(x_ss(2),'.')
subplot(3,2,3)
boxplot(xAll(3,:), 'TLC Wireless')
hold on
plot(x_ss(3),'.')
subplot(3,2,4)
boxplot(xAll(4,:), 'TLC Wired')
hold on
plot(x_ss(4),'.')
subplot(3,2,5)
boxplot(xAll(5,:), 'Water Management')
hold on
plot(x_ss(5),'.')
subplot(3,2,6)
boxplot(xAll(6,:), 'Rail Transportation')
hold on
plot(x_ss(6),'.')

figure (10)
subplot(3,2,1)
boxplot(xAll(7,:), 'Finance')
hold on
plot(x_ss(7),'.')
subplot(3,2,2)
boxplot(xAll(8,:), 'Fuel&Petro Grid')
hold on
plot(x_ss(8),'.')
subplot(3,2,3)
boxplot(xAll(9,:), 'Natural Gas')
hold on
plot(x_ss(9),'.')
subplot(3,2,4)
boxplot(xAll(10,:), 'Naval Ports')
hold on
plot(x_ss(10),'.')
subplot(3,2,5)
boxplot(xAll(11,:), 'Satellite Com&Navigation')
hold on
plot(x_ss(11),'.')
    
% ======================
% Question 2.17 - Uncertainty in External Shocks
% ======================
% Assume values of A are now random variables
% Uniform distribution +/- 10% of give values
% Use Monte Carlo to determine uncertainty in x given uncertainty in A

% HEM 11/10: deleted old code

runs = 10000; % number of MC runs
% Initialize storage matrix for f values
xAll_randf = NaN(m,runs);
% Use snowstorm external shock f_ss 
for i = 1:runs
    % Generate vector of random multipliers
    % Want random values between 0.9 and 1.1, inclusive
    mults = (1.1-0.9) .* rand(m,1) + 0.9;
    % Use element-wise vector multiplication to sample a new vector f
    f_samp = f_ss .* mults; 
    % Calculate x vector 
    x_samp = (I-A)\f_samp;
    % Put into the storage matrix
    xAll_randf(:,i) = x_samp;
end

% HEM 11/10 : added documentation changed code to match new
% plot the distributions of uncertainty for each sector in histogram
% with a fit curve to show the probability density function
figure(11)
%suptitle ('Uncertainty distributions of degradation for sectors 1-6')
subplot(3,2,1)
histfit(xAll_randf(1,:));
title('Air Transportation');
subplot(3,2,2)
histfit(xAll_randf(2,:));
title('Electricity');
subplot(3,2,3)
histfit(xAll_randf(3,:));
title('Wireless Telecom');
subplot(3,2,4)
histfit(xAll_randf(4,:));
title('Wired Telecom');
subplot(3,2,5)
histfit(xAll_randf(5,:));
title('Water Management');
subplot(3,2,6)
histfit(xAll_randf(6,:));
title('Rail Transportation');

figure (12)
%suptitle('Uncertainty distributions for sectors 7-11')
subplot (3,2,1)
histfit(xAll_randf(7,:));
title('Finance');
subplot(3,2,2)
histfit(xAll_randf(8,:));
title('Fuel Grid');
subplot(3,2,3)
histfit(xAll_randf(9,:));
title('Natural Gas')
subplot(3,2,4)
histfit(xAll_randf(10,:));
title('Naval Ports')
subplot(3,2,5)
histfit(xAll_randf(11,:));
title('Satellite Com')

%HEM 11/10 : changed variable names and documentation
%calculate the mean values
xMean_randf = mean(xAll_randf,2);  
diff_randf = xMean_randf - x_ss  ;
%compare the difference for the two uncertainty analyses
disp(sum(abs(diff)))
disp(sum(abs(diff_randf)))

%calculate the standard deviation
stdev_randf = std(xAll_randf')'; 
%compare the standard deviations for the two uncertainty analyses
disp(sum(stdev))
disp(sum(stdev_randf))

%Use a box plot to represent uncertainty
figure (13)
subplot(3,2,1)
boxplot(xAll_randf(1,:), 'Air Transportation')
hold on
plot(x_ss(1),'.')
subplot(3,2,2)
boxplot(xAll_randf(2,:), 'Electricity')
hold on
plot(x_ss(2),'.')
subplot(3,2,3)
boxplot(xAll_randf(3,:), 'TLC Wireless')
hold on
plot(x_ss(3),'.')
subplot(3,2,4)
boxplot(xAll_randf(4,:), 'TLC Wired')
hold on
plot(x_ss(4),'.')
subplot(3,2,5)
boxplot(xAll_randf(5,:), 'Water Management')
hold on
plot(x_ss(5),'.')
subplot(3,2,6)
boxplot(xAll_randf(6,:), 'Rail Transportation')
hold on
plot(x_ss(6),'.')

figure (14)
subplot(3,2,1)
boxplot(xAll_randf(7,:), 'Finance')
hold on
plot(x_ss(7),'.')
subplot(3,2,2)
boxplot(xAll_randf(8,:), 'Fuel&Petro Grid')
hold on
plot(x_ss(8),'.')
subplot(3,2,3)
boxplot(xAll_randf(9,:), 'Natural Gas')
hold on
plot(x_ss(9),'.')
subplot(3,2,4)
boxplot(xAll_randf(10,:), 'Naval Ports')
hold on
plot(x_ss(10),'.')
subplot(3,2,5)
boxplot(xAll_randf(11,:), 'Satellite Com&Navigation')
hold on
plot(x_ss(11),'.')

% ======================
% Question 2.18 - Comparing Lower and Upper Bounds
% ======================
% Compare results from running lower and upper bounds for both A and f_ss

% lower bound for A (-10%)
A_lower = 0.9 .* A; % SC 11/06: dot multiply for correctness
x_lowerA = (I-A_lower)\f_ss; % SC 11/06: change to mldivide

% lower bound for f (-10%)
f_lower = 0.9 .* f_ss;
x_lowerf = (I-A)\f_lower;

% upper bound for A (+10%)
A_upper = 1.1 .* A;
x_upperA = (I-A_upper)\f_ss;

% upper bound for f (+10%)
f_upper = 1.1 .* f_ss;
x_upperf = (I-A)\f_upper;


% SC 11/06: didn't touch below
    %make comparison
figure(15)
x=1:1:11;
plot(x,x_lowerA,x,x_upperA)
xlabel('sector')
ylabel('degradation')
title('Degradation with the lower and the upper bound of uncertainty in A')
legend('Lower Bound','Upper Bound')

figure(16)
x=1:1:11;
plot(x,x_lowerf,x,x_upperf)
xlabel('sector')
ylabel('degradation')
title('Degradation with the lower and the upper bound of uncertainty in f')
legend('Lower Bound','Upper Bound')

% HEM 11/10: updated code and documentation
% ======================
% Question 2.19 - maximum impacts absorption
% ======================

% Maximize f subjected to these constraints
% degradation to fuel petroleum grid, TLC wireless and finance cannot
% be greater than 20%, 30% and 35%:
% S(8,:)*f<=0.2
% S(3,:)*f<=0.3
% S(7,:)*f<=0.35
% Each element in f should be positive: -fi<=0
% Optimization method using the linprog function
% Minimize f_opt = - f 

f_opt=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]';

options = optimoptions('linprog', 'Algorithm', 'dual-simplex');

A_opt=[S(8,:); S(3,:); S(7,:)];
B_opt=[0.2;0.3;0.35];
Aeq = [];
beq = [];
lb= zeros(11,1); %lower bound is zero
ub= ones(11,1); %upper bound is one
                          
% Run different versions depending on Matlab version
verStr = version;
verS = strtok(verStr);
ver = str2double(verS(1:3));

if ver >= 9.1 % version 2016
    [x, fval, exitflag, output, lambda] = ...
        linprog(f_opt, A_opt, B_opt, Aeq, beq, lb, ub, options);
else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, fval, exitflag, output, lambda] = ...
        linprog(f_opt, A_opt, B_opt, Aeq, beq, lb, ub, [], options);
    % Basically adds a x0 vector input before options
end