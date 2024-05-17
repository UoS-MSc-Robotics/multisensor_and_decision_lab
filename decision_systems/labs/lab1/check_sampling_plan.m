% Analyze different sampling plans
clc; clear; close all;

% Add sampling and evaluation functions to the path
addpath(fullfile(pwd, '/lab1/sampling/'));
addpath(fullfile(pwd, '/lab1/evaluation/'));

% Labels for design constraints
design_constraints = {'max pole', 'gain margin', 'phase margin', 'rise time', 'peak time', 'overshoot', ...
    'undershoot', 'settling time', 'steady-state error', 'control input'};

% Sampling plans to analyze
sampling_plans_list = {'full factorial', 'sobol set', 'latin hypercube', 'random Latin hypercube'};

% Analyze different sampling plans
[P, best_sampling_plan] = mmphi_analysis(sampling_plans_list);
fprintf('The best sampling plan is: %s\n', best_sampling_plan);

% Evaluate the control system
Z = evaluateControlSystem(P);

% Specify the goals
goals = [1, 6, 70, 2, 10, 10, 8, 20, 1, 0.67];

% Implement knowledge discovery
knowledge_discovery(Z, best_sampling_plan, design_constraints, goals);



% Function to analyze different sampling plans using mmphi
function [P_best, best_sampling_plan] = mmphi_analysis(sampling_plans_list)
    scale = 1;
    q = [10, 10];
    Edges = 1;
    min_phi_metric = Inf;
    n_rows = ceil(length(sampling_plans_list) / 2);

    figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    for i = 1:length(sampling_plans_list)
        sampling_plan = sampling_plans_list{i};
        if strcmp(sampling_plan, 'full factorial')
            P = fullfactorial(q, Edges);
        elseif strcmp(sampling_plan, 'sobol set')
            P = sobolset(length(q));
            P = net(P, q(1)*q(2));
        elseif strcmp(sampling_plan, 'latin hypercube')
            P = lhsdesign(q(1)*q(2), length(q));
        elseif strcmp(sampling_plan, 'random Latin hypercube')
            P = rlh(q(1)*q(2), length(q), Edges);
        else
            error('Invalid sampling plan specified.');
        end

        phi_metric = mmphi(P * scale, 2, 2);
        fprintf('The MMPhi metric for %s sampling plan is: %f\n', sampling_plan, phi_metric);

        if phi_metric < min_phi_metric
            min_phi_metric = phi_metric;
            best_sampling_plan = sampling_plan;
            P_best = P;
        end

        % Plot the sampling plan using subplot
        subplot(n_rows, 2, i);
        plot(P(:,1), P(:,2), 'o');
        title(sprintf('%s', sampling_plan));
        xlabel('x_1');
        ylabel('x_2');
    end
end


% Function to implement knowledge discovery
function knowledge_discovery(Z, best_sampling_plan, design_constraints, goals)
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    stability1 = nan * ones(size(Z, 2), 1);
    stability2 = nan * ones(size(Z, 2), 1);

    stability1(1:3) = goals(1:3);
    stability2(1:3) = goals(1:3);
    stability2(3) = 30;

    transient = nan * ones(size(Z, 2), 1);
    transient(4:8) = goals(4:8);

    steady_state = nan * ones(size(Z, 2), 1);
    steady_state(9) = goals(9);

    sustainability = nan * ones(size(Z, 2), 1);
    sustainability(10) = goals(10);


    % Append extra line to A
    Z = [Z; stability1'; stability2'; transient'; steady_state'; sustainability'];

    % Grouping vector
    % groupDataVector = ones(1, size(Z, 1) - 2);
    % groupDataVector = [groupDataVector, 2, 2];

    % create a vector of strings
    groupDataVector = cell(1, size(Z, 1) - 5);
    for i = 1:length(groupDataVector)
        groupDataVector{i} = 'All Design Evaluation';
    end

    groupDataVector = [groupDataVector, 'Stability', 'Stability', 'Transient', 'Steady State', 'Sustainability'];

    p = parallelplot(Z, 'groupData', groupDataVector, 'Color', {'green','red', 'blue', 'magenta', 'black'}, 'LineWidth', 2);
    p.MarkerStyle = {'none','o', 'o', 'o', 'o'};
    p.MarkerSize(end) = 15;
    p.LineStyle = {'-', '--', '--', '--', '--'};
    p.LegendTitle = 'Performance Criteria';


    p.CoordinateTickLabels = design_constraints;
    p.YLabel = 'Performance metric value';
    p.Title = sprintf('Performance evaluations for %s sampling plan', best_sampling_plan);


end
