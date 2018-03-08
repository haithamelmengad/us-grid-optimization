% S Cai, H El Mengad, K Xie
% CEE 498 SIS Project Part 4

% This file contains the problem data for Part 4.38
% Run this code to save the data to memory.
% Everything is a column vector for now.

% -------------------
% Table 3: Generating Plant Production Data
% Plant numbers for reference
% e = existing
% n = new, can build to selected capacity
% 1 e conventional coal
% 2 e advanced coal
% 3 e conventional combustion nat gas
% 4 e nuclear
% 5 n advanced coal with ccs
% 6 n advanced CC2 nat gas
% 7 n wind
% 8 n solar
% 9 e hydroelectric

% xMax_i = maximum generating capacity of plant i (MW)
xMax_i = [250, 210, 180, 400, 500, 500, 300, 400, 200]';

% xMin_i = minimum (must run) capacity of plant i (MW)
xMin_i = [70, 0, 0, 100, 0, 0, 0, 0, 0]';

% ou_i = unplanned outage rate
ou_i = [0.05, 0.04, 0.035, 0.05, 0.05, 0.05, 0.03, 0.03, 0.08]';

% op_i = planned osdutage rate
op_i = [0.15, 0.15, 0.70, 0.10, 0.15, 0.15, 0.66, 0.80, 0.12]';

% g_i = global warming potential of plant i (kg CO2e/MWh)
g_i = [1001, 766, 443, 16, 16, 305, 12, 45, 78]';


% -------------------
% Table 4: Generating Plant Cost Data

% civBar = avg variable operating cost for plant i ($/MWh)
civBar = [28.6, 29.1, 79.9, 11.5, 36.8, 44.3, 0.9, 1, 1.5]';

% civSD = SD for variable op cost for plant i ($/MWh)
civSD = [3, 3, 3, 4, 2, 2, 0.15, 0.5, 0.2]';

% cicBar = avg annualized capital op cost for plant i ($/kW)
% NaN for existing plants
cicBar = [NaN, NaN, NaN, NaN, 43.3, 17.8, 83.3, 145, NaN]';

% cicSD = SD for annualized capital op cost for plant i ($/kW)
% NaN for existing plants
cicSD = [NaN, NaN, NaN, NaN, 9, 6, 12, 20, NaN]';


% -------------------
% Table 5: Load Blocks for annual power demand
% All 8766 hours in a year divided into 6 load blocks
% Load blocks indexed by 't'

% n_t = number of hours in load block t (hr)
n_t = [100, 315, 671, 1245, 2780, 3655]';

% l_t = load in load block t (MW)
l_t = [1390, 1305, 830, 545, 355, 200]';


% -------------------
% Table 6: Demand side management programs
% 3 DSM programs, indexed by 'k'

% sMax_kt = max energy savngs by program k in block t (MW)
% k by t matrix (3x6)
% 1 energy efficiency
% 2 energy efficiency
% 3 load control
sMax_kt = [80, 30, 25, 12, 5, 0;
    70, 30, 15, 10.5, 5, 0;
    100, 0, 0, 0, 0, 0];

% ckdBar = avg cost of DSM program k implementation ($/MWh)
ckdBar = [55, 65, 100]';

% ckdSD = cost SD of DSM program k implementation ($/MWh)
ckdSD = [15, 7, 20]';

% Decision variable sets
I = 9;
T = 6;
K = 3;

% % ---------Constraints-----------
% ---------------------------------
% LINEAR INEQUALITIES Ax <= b
% ---------------------
% ---------------------
% Total 161 constraints consisting of:
% i Load constraints (6)
% ii Instantaneous capacity constraints (54)
% iii Annual energy constraints existing plants (5)
% iv Annual energy constraints new plants (4)
% v Minimum generation constraints (54)
% vi Bounds on DSM programs (6)
% vii New generation bounds (32)

% However, ii, v, vi, and last 8 of vii are lower & upper bounds
% See below in bounds section for these formulations
% 39 constraints in this section
% 122 constraints in UB/LB section

% Initialize a 39x61 matrix for A
% Initialize a 39x1 vector for b
A = NaN(39, 61);
b = NaN(39, 1);

% ---------------------
% i Load constraints (6)
% l_t - sum_k z_k sMax_kt <= sum_i x_it, for all T
% one constraint for each value of t
% x_it = -1 for all i for each t
% i.e. every 6 elements in a row are -1

for t = 1:T % first 6 rows
    for i = 1:54 % set remaining x_it
        if mod(i-t, 6) == 0
            A(t,i) = -1;
        else
            A(t,i) = 0;
        end
    end
end

% All y_i constants are zero
A(1:6, 55:58) = 0;

% All z_k constants are -1* sumK sMax_kt for each value of t
for t = 1:T
    ctr = 59;
    for k = 1:K
        A(t, ctr) = -1 * sMax_kt(k,t);
        ctr = ctr + 1;
    end
    b(t) = -1 * l_t(t);
end

% ---------------------
% iii annual energy constraints existing (5)
% sumT n_t x_it <= (1-op_i) 8766 xMax_i
% one line for each i = 1, 2, 3, 4, 9 
existing_i = [1, 2, 3, 4, 9];

% Initialize all constants to be zero
% Includes all y_i, z_k, and non-relevant values of x_it
A(7:11,:) = 0;

% set x_it values
for ctr = 7:11 % next 5 contraints
    % set value of i
    i = existing_i(ctr - 6); % 5 values only
    % calculate the indices of x correspond to input i 
    inds = (i*6-5):(i*6);
    for t = 1:T
        A(ctr, inds(t)) = n_t(t);
    end
    % set values of b for each i
    b(ctr) = (1 - op_i(i)) * 8766 * xMax_i(i);
end

% ---------------------
% iv annual energy constraints new (4)
% sumT n_t x_it - (1-op_i)8766 * y_i <= 0
% i = 5, 6, 7, 8
% rows 12 to 15 in matrix A
new_i = 5:8;

A(12:15,:) = 0;


for ctr = 12:15 % next 5 contraints
    % set value of i
    i = new_i(ctr - 11); % 4 values only
    % set x_it values
    % calculate the indices of x correspond to input i 
    inds = (i*6-5):(i*6);
    for t = 1:T
        A(ctr, inds(t)) = n_t(t);
    end
    % set the single y_i value
    % y_i's are 55-58
    A(ctr, 43+ctr) = -8766 * (1 - op_i(i));
    % set values of b
    b(ctr) = 0; 
end

% ---------------------
% vii new generation bounds (first 24)
% rows 16-39 in matrix A
% x_it - (1-ou_i)*y_i <= 0
% for all t
% i in [5,8]
% single x_it and single y_i value in each row
% z_k is zero

A(16:39,:) = 0;
ctr = 16;
% find columns for x_i and y_i
ycol = 55; 
xcol = 25; % first 4 * 6 = 25 are zero

for i = 5:8
    for t = 1:T
        A(ctr, xcol) = 1;
        A(ctr, ycol) = -1 * (1-ou_i(i));
        xcol = xcol + 1;
        ctr = ctr + 1;
    end
    ycol = ycol + 1;
end

% values of b are zero
b(16:end) = 0;

% ---------------------
% ---------------------
% LINEAR EQUALITIES Aeq*x <= beq
% ---------------------
% ---------------------
% There are no linear equalities
Aeq = [];
beq = [];


% ---------------------
% ---------------------
% BOUNDS lb <= x <= ub
% ---------------------
% ---------------------
% Constraints that are bounds on decision variables
% ii Instantaneous capacity constraints (54)
% v Minimum generation constraints (54)
% vi Bounds on DSM programs (6)
% vii New generation bounds (last 8)

lb = NaN(61, 1);
ub = NaN(61, 1);

% Bounds of x_it (first 54)
% ii Instantaneous capacity constraints
% x_it <= (1-ou_i) * xMax_i
% v Minimum generation constraints
% x_it >= x_i^min
% x11 = x1min, x12 = x1min, etc. so blocks of 9 for both
% IE lower bound for 

ctr = 1;
for i = 1:I
    for t = 1:T
        lb(ctr) = xMin_i(i);
        ub(ctr) = (1 - ou_i(i)) * xMax_i(i);
        ctr = ctr + 1;
    end
end


% Bounds on y_i (next 4)
% vii new generation bounds for i in [5,8]
% y_i >= 0
% y_i <= xMax_i

for i = 5:8
    lb(ctr) = 0;
    ub(ctr) = xMax_i(i);
    ctr = ctr + 1;
end


% Bounds on z_k (last 3)
% vi bounds on dsm programs
% z_k >= 0
% z_k <= 1

for k = 1:K
    lb(ctr) = 0;
    ub(ctr) = 1;
    ctr = ctr + 1;
end
% --------Constraints End-----------
% ----------------------------------