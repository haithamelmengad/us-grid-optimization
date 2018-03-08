% SC, HEM, KX
% 498 SIS Project Part 4.32
% Add Constraints of tax on emissions ($/kg CO2eq)

% tax on CO2 eq emissions is being considered 
% in the minimum cost optimization problem
% Find: under what values will the emissions tax influence the solution?
% New objective function: f = f1 (the cost function in pt428) 
%              + 1000 * tax * f2 (the emission function in 429)

% first load the problem setup data
load all_problem_data

%% get f1 from 428
f = NaN(61,1);
ctr = 1;
for i = 1:I
    for t = 1:T
        f(ctr) = civBar(i) * n_t(t);
        ctr = ctr + 1;
    end
end

for i = 5:8
    f(ctr) = 1000 * cicBar(i);
    ctr = ctr + 1;
end

for k = 1:K
    f(ctr) = ckdBar(k) * (sMax_kt(k,:) * n_t);
    ctr = ctr + 1;
end
f1 = f;

% get f2 from 429
f = NaN(61,1);

ctr = 1;
for i = 1:I
    for t = 1:T
        f(ctr) = 0.001 * g_i(i) * n_t(t);
        ctr = ctr + 1;
    end
end

for i = 5:8
    f(ctr) = 0;
    ctr = ctr + 1;
end

for k = 1:K
    f(ctr) = 0;
    ctr = ctr + 1;
end
f2 = f;

%% Simulate different tax values
N = 0.01; 
numSteps = length(0:N:20);

% Initialize a cell array to store linprog results
TaxRes = cell(5, numSteps);
% Storage matrix for tax and min cost results
fTax = NaN(numSteps,3);
count = 1;

for tax = 0:N:20 % simulate tax values from 0~$20
    
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex');

% Run different versions depending on Matlab version
verStr = version;
verS = strtok(verStr);
ver = str2double(verS(1:3));

if ver >= 9.1 % version 2016
    [x, fval, exitflag, output, lambda] = ...
        linprog(f1+1000*tax*f2, A, b, Aeq, beq, lb, ub, options);
else % for Ketong's 2014. assuming anything before 2016 is like this
    [x, fval, exitflag, output, lambda] = ...
        linprog(f1+1000*tax*f2, A, b, Aeq, beq, lb, ub, [], options);
    % Basically adds a x0 vector input before options
end
     % Save values
    TaxRes{i, 1} = x;
    TaxRes{i, 2} = fval;
    TaxRes{i, 3} = exitflag;
    TaxRes{i, 4} = output;
    TaxRes{i, 5} = lambda;
    
    % Store both tax and min cost in a matrix for charting
    fTax(count,:) = [tax, fval, f2'*x];
    count = count + 1;
end


% Horizontal Lines at the Min Cost and Min GWP
minCostLine = zeros(size(fTax,1), 1);
minCostLine(:) = 5.5031E07;

minGWPLine = zeros(size(fTax,1), 1);
minGWPLine(:) = 6.6283E05;

% plot the results, cap and the fval
figure(1)
plot(fTax(:,1),fTax(:,2),'.',fTax(:,1),minCostLine,'-')
% axis([0 2 4e7 15e8]);
xlabel('Tax Value')
ylabel('Expected Cost')
title('Tax Value Embedded in Cost Optimization Problem')

% OK, tax of $20 was too high, zoom in to sub 50c

figure(9)
plot(fTax(:,1),fTax(:,2),'.',fTax(:,1),...
    minCostLine,'-','markers',20, 'LineWidth', 2)
axis([0 0.5 0 5e8]);
xlabel('Tax Value ($/kg CO2e)', ...
    'FontSize', 14)
ylabel('Expected Cost ($)', ...
    'FontSize', 14)
hl = legend('Cost', 'Minimum Cost');
set(hl, 'FontSize', 14);
title('Tax Value Embedded in Cost Optimization Problem', ...
    'FontSize', 16, 'fontWeight', 'bold')

% the result shows that only when the tax is set as 0 
% will the minimum cost optimization problem not be influenced


%% Find: choose a tax which will result in the same solution as the minimum 
% emissions optimization problem
figure(3)
plot(fTax(:,1),fTax(:,3),'.',...
    fTax(:,1), minGWPLine, '-', 'markers', 10)
axis([0 12 6e5 9e5])
xlabel('Tax Value ($/kg CO2e)')
ylabel('Emissions (MT CO2e)')
legend('Emissions', 'Minimum Emissions')
title('Tax Value Embedded in Emission Optimization Problem')

% the result shows that the tax should be set as $10.91 in order to
% meet the same solution of the emission opt problem
