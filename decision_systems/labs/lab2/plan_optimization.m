% Analyze different sampling plans
clc; clear; close all;

% Add sampling and evaluation functions to the path
addpath(fullfile(pwd, '/lab1/sampling/'));
addpath(fullfile(pwd, '/lab1/evaluation/'));

% Sampling plans to analyze
sampling_plans_list = {'full factorial', 'sobol set', 'latin hypercube', 'random Latin hypercube'};

% Analyze different sampling plans
[P, best_sampling_plan] = mmphi_analysis(sampling_plans_list);
fprintf('The best sampling plan is: %s\n', best_sampling_plan);

% Optimize the sampling plan
Z = optimizeControlPlan(P);

% Implement knowledge discovery
knowledge_discovery(Z, best_sampling_plan);



% Function to optimize the sampling plan
function Z_optimized = optimizeControlPlan(P)
    Z = evaluateControlSystem(P);

    % Step 1: Convert gain margin to decibels
    Z(:,2) = 20*log10(Z(:,2));

    % Step 2: Minimize the gain margin
    Z(:,2) = -Z(:,2);

    % Step 3: Minimize the absolute deviation of the gain margin within 30 and 70 dB
    Z(:,3) = abs(Z(:,3) - 50);

    % Step 4: Remove all inf values in Z, set to a high value
    Z(isinf(Z)) = 1e3;

    % Return the optimized sampling plan
    Z_optimized = Z;

end


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
            X = fullfactorial(q, Edges);
        elseif strcmp(sampling_plan, 'sobol set')
            X = sobolset(length(q));
            X = net(X, q(1)*q(2));
        elseif strcmp(sampling_plan, 'latin hypercube')
            X = lhsdesign(q(1)*q(2), length(q));
        elseif strcmp(sampling_plan, 'random Latin hypercube')
            X = rlh(q(1)*q(2), length(q), Edges);
        else
            error('Invalid sampling plan specified.');
        end

        phi_metric = mmphi(X * scale, 5, 1);
        fprintf('The MMPhi metric for %s sampling plan is: %f\n', sampling_plan, phi_metric);

        if phi_metric < min_phi_metric
            min_phi_metric = phi_metric;
            best_sampling_plan = sampling_plan;
            P_best = X;
        end

        % Plot the sampling plan using subplot
        subplot(n_rows, 2, i);
        plot(X(:,1), X(:,2), 'o');
        title(sprintf('%s', sampling_plan));
        xlabel('x_1');
        ylabel('x_2');
    end
end



% Function to implement knowledge discovery
function knowledge_discovery(Z, best_sampling_plan)
    labels = {'max pole', 'gain margin', 'phase margin', 'rise time', 'peak time', 'overshoot', 'undershoot', 'settling time', 'steady-state error', 'control input'};

    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    p = parallelplot(Z, 'Color', 'b');
    p.CoordinateTickLabels = labels;
    % set y-axis labels
    p.YLabel = 'Performance metric value';
    p.Title = sprintf('Performance evaluations for %s sampling plan', best_sampling_plan);
end
