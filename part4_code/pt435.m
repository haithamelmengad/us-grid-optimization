% pt435.m
% Proposed Electric Utility Design
% CEE 498 SIS Project
% S Cai, K Xie, H El Mengad
% Contains the triple-objective optimization only.

% ---------------------------
% Summary of prior results
% 4.28 Min Cost: $5.5E07, emissions 1.03E06 MT
% 4.29 Min Emissions: $1.78E08, 6.63E05 MT
% 4.30 Cost vs. Emissions optimals: several
% inflection points, bottom left-most was
% $8E07 and <7E05 MT
% 4.31 Emissions cap: couple inflection points,
% one matched $8E07 and 7F05 MT from 4.30
% 4.32 Emissions tax: $1.2/kg basically caused
% min emissions. If we select target emissions of
% 7E5 we will be prepared for carbon tax
% around 20c per kg, which is order of mag
% above many current carbon taxes
% 4.33 Min Var: variance minimized and relatively
% flat change in cost range $7.582-7.588E07,
% so selection of $8E07 cost is fine for min var

% So based on prior analyses, 
% goal is to target ~$8E07 and 7E05 MT.
% However, want to check this by completing
% some additional optimizations

% ---------------------------
% Structure:
% load the problem setup data
load all_problem_data

% 1. Perform triple-objective opt. btwn
% cost, GHG, and variance to make sure

%% DEFINE OBJECTIVE FUNCTIONS FOR TRIPLE OBJECTIVE OPTIMIZATION
% Use quadprog, so need to set up H matrices for all 3 problems
% Use two factors, theta and phi, third is 1-theta-phi

% ---------------------------
% ---------------------------
% i. Set up min expected cost quadprog
% Use 4.34 code

% get H1 and f1 for cost, see pt428
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



% ---------------------------
% ---------------------------
% ii. Set up min GHG quadprog
% Use 4.29 code

% H matrix is zeros again
H2 = zeros(61);

% f setup vector is taken from 4.29
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


% ---------------------------
% ---------------------------
% iii. Set up min variance quadprog
% Use 4.34 code

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

% This is the full matrix going into the quadprog
% which is HT E H
Hvar = 2 * H3' * sigma * H3;



%% SETUP FOR TRIPLE OBJECTIVE OPT

% Set normalization factors based on results of 4.28, 4.29, 4.33
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

% Initialize a cell array to store quadprog results
tripleObjRes = cell(numRuns, 5);

% Storage matrix for f1' and f2' results
% Cols 1-3 are to list our values of 3 factors
% Cols 4-6 are to store the fvals
fvalComp = NaN(numRuns,5);


%% RUN TRIPLE OBJECTIVE OPTIMIZATION

    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');

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
        quadprog(HTemp, fTemp, A, b, Aeq, beq, lb, ub, [], options);

    % Save values
    tripleObjRes{i, 1} = x;
    tripleObjRes{i, 2} = fval;
    tripleObjRes{i, 3} = exitflag;
    tripleObjRes{i, 4} = output;
    tripleObjRes{i, 5} = lambda;
    
    % Store results in charting matrix
    fvalComp(i,1) = fac1;
    fvalComp(i,2) = fac2;
    fvalComp(i,3) = fac3;
    % These are individual calcs
    fvalComp(i,4) = (0.5*x'*H1*x) + (f1'*x);
    fvalComp(i,5) = (0.5*x'*H2*x) + (f2'*x);
    fvalComp(i,6) = (0.5*x'*Hvar*x) + (f3'*x);
end



%% Plot the results

% subplots
figure(5)
subplot(2,2,1)
% 3D visualization
scatter3(fvalComp(:, 4), fvalComp(:, 5), fvalComp(:, 6))
rotate3d on
xlabel('Min Expected Cost')
ylabel('Min GHG')
zlabel('Min Cost Variance')
title('Triple-Objective Optimization')

% 2D of cost-ghg with var color
subplot(2,2,2)
scatter(fvalComp(:, 4), fvalComp(:, 5), [], fvalComp(:, 6))
colorbar
xlabel('Min Expected Cost')
ylabel('Min GHG')
ylabel(colorbar, 'Variance')
title('Cost vs. GHG with Variance colors')

% 2D of cost-var with ghg color
subplot(2,2,3)
scatter(fvalComp(:, 4), fvalComp(:, 6), [], fvalComp(:, 5))
colorbar
xlabel('Min Expected Cost')
ylabel('Min Var')
ylabel(colorbar, 'GHG')
title('Cost vs. Var with GHG colors')

% 2D of ghg-var with cost color
subplot(2,2,4)
scatter(fvalComp(:, 5), fvalComp(:, 6), [], fvalComp(:, 4))
colorbar
xlabel('Min GHG')
ylabel('Min Var')
ylabel(colorbar, 'Expected Cost')
title('GHG vs. Var with Cost colors')


% individual plots
figure(6)
% 3D visualization
scatter3(fvalComp(:, 4), fvalComp(:, 5), fvalComp(:, 6))
rotate3d on
xlabel('Min Expected Cost')
ylabel('Min GHG')
zlabel('Min Cost Variance')
title('Triple-Objective Optimization')

% 2D of cost-ghg with var color
figure(7)
scatter(fvalComp(:, 4), fvalComp(:, 5), [], fvalComp(:, 6))
colorbar
xlabel('Min Expected Cost')
ylabel('Min GHG')
ylabel(colorbar, 'Variance')
title('Cost vs. GHG with Variance colors')

% 2D of cost-var with ghg color
figure(8)
scatter(fvalComp(:, 4), fvalComp(:, 6), [], fvalComp(:, 5))
colorbar
xlabel('Min Expected Cost')
ylabel('Min Var')
ylabel(colorbar, 'GHG')
title('Cost vs. Var with GHG colors')

% 2D of ghg-var with cost color
figure(9)
scatter(fvalComp(:, 5), fvalComp(:, 6), [], fvalComp(:, 4))
colorbar
xlabel('Min GHG')
ylabel('Min Var')
ylabel(colorbar, 'Expected Cost')
title('GHG vs. Var with Cost colors')


