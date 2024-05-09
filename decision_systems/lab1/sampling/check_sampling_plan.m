% Analyze different sampling plans
close all;

scale = 1;

% Full factorial sampling plan
% Generates a full factorial sampling plan in the unit cube
%
% Inputs:
%       q - k-vector containing the number of points along each dimension
%       Edges - if Edges=1 the points will be equally spaced from edge to
%               edge (default), otherwise they will be in the centres of 
%               n = q(1)*q(2)*...q(k) bins filling the unit cube.
%
% X - full factorial sampling plan

q = [10 10];
Edges = 1;
X_ff = fullfactorial(q, Edges) * scale;

% Plot the full factorial sampling plan
figure;
plot(X_ff(:,1), X_ff(:,2), 'o');
title('Full factorial sampling plan');
xlabel('x_1');
ylabel('x_2');


% Random Latin hypercube sampling plan
% Generates a random Latin hypercube within the [0,1]^k hypercube.
%
% Inputs:
%       n - desired number of points
%       k - number of design variables (dimensions)
%       Edges - if Edges=1 the extreme bins will have their centres on the
%               edges of the domain, otherwise the bins will be entirely 
%               contained within the domain (default setting). 
%
% Output:
%       X - Latin hypercube sampling plan of n points in k dimensions.
q = [10 10];
Edges = 1;
X_rlh = rlh(q(1)*q(2), length(q), Edges) * scale;

% Plot the random Latin hypercube sampling plan
figure;
plot(X_rlh(:,1), X_rlh(:,2), 'o');
title('Random Latin hypercube sampling plan');
xlabel('x_1');
ylabel('x_2');

% Analyze the sampling plan using mmphi

phi_metric_ff = mmphi(X_ff, 5, 1);
phi_metric_rlh = mmphi(X_rlh, 5, 1);

disp(['The value of the phi metric for the full factorial sampling plan is ', num2str(phi_metric_ff)]);
disp(['The value of the phi metric for the random Latin hypercube sampling plan is ', num2str(phi_metric_rlh)]);


% Knowledge discovery
%
% Simulation model for a PI controller.
% ACS6124 Part 2 2023/24
% RP, 6 February 2024
%
% Z = evaluateControlSystem(X);
%
% Input:  X - a sampling plan
%             (one design per row, one design variable per column)
% Output: Z - performance evaluations
%             (one design per row, one criterion per column;
%              criteria are...
%              1:  maximum closed-loop pole magnitude
%              2:  gain margin
%              3:  phase margin
%              4:  10-90% rise time
%              5.  peak time
%              6.  overshoot (% points)
%              7.  undershoot (% points)
%              8.  2% settling time
%              9.  steady-state error (% points))
%              10. aggregate control input (MJ)
%

Z_ff = evaluateControlSystem(X_ff);
Z_rlh = evaluateControlSystem(X_rlh);

labels = {'max pole', 'gain margin', 'phase margin', 'rise time', 'peak time', 'overshoot', 'undershoot', 'settling time', 'steady-state error', 'control input'};

figure;
p = parallelplot(Z_ff, 'Color', 'b');
p.CoordinateTickLabels = labels;
% set y-axis labels
p.YLabel = 'Performance metric value';
p.Title = 'Performance evaluations for full factorial sampling plan';

figure;
p = parallelplot(Z_rlh, 'Color', 'r');
p.CoordinateTickLabels = labels;
% set y-axis labels
p.YLabel = 'Performance metric value';
p.Title = 'Performance evaluations for random Latin hypercube sampling plan';