% Copyright (c) 2024 Leander Stephen D'Souza

% Program to analyze different sampling plans and build an optimizing engine

clc; clear; close all;

% Build all the mex files
mex(fullfile(pwd, '/EA_Toolbox/rank_nds.c'));
mex(fullfile(pwd, '/EA_Toolbox/crowdingNSGA_II.c'));
mex(fullfile(pwd, '/EA_Toolbox/btwr.c'));
mex(fullfile(pwd, '/EA_Toolbox/sbx.c'));
mex(fullfile(pwd, '/EA_Toolbox/polymut.c'));
mex('-DVARIANT=4', ...
    fullfile(pwd, '/Hypervolume/Hypervolume_MEX.c'), ...
    fullfile(pwd, '/Hypervolume/hv.c'), ...
    fullfile(pwd, '/Hypervolume/avl.c'));
mex(fullfile(pwd, '/EA_Toolbox/rank_prf.c'));

% Add sampling and evaluation functions to the path
addpath(fullfile(pwd, '/sampling/'));
addpath(fullfile(pwd, '/evaluation/'));
addpath(fullfile(pwd, '/EA_Toolbox/'));

% Labels for the design constraints
design_constraints = {'max pole', 'gain margin', 'phase margin', 'rise time', 'peak time', 'overshoot', ...
    'undershoot', 'settling time', 'steady-state error', 'control input'};

% Sampling plans to analyze
sampling_plans_list = {'full factorial', 'sobol set', 'latin hypercube', 'random Latin hypercube'};

% Analyze different sampling plans
[P, best_sampling_plan] = mmphi_analysis(sampling_plans_list);
fprintf('The best sampling plan is: %s\n', best_sampling_plan);

% % Find the lower dimensional representation of the data
data_mining(evaluateControlSystem(P));

% Building the optimizing engine
% Initial Client Assignment
iterations = 50;
priority = [3, 2, 2, 1, 0, 1, 0, 0, 1, 2];
weighted_priority = [0.8, 0.6, 0.6, 0.4, 0.2, 0.4, 0.2, 0.2, 0.6, 0.6];
goals = [1, -6, 20, 2, 10, 10, 8, 20, 1, 0.67];

buildOptimizingEngine(true, P, iterations, goals, priority, weighted_priority, best_sampling_plan, design_constraints);

% Reduce control input
goals(10) = 0.63;
buildOptimizingEngine(true, P, iterations, goals, priority, weighted_priority, best_sampling_plan, design_constraints);



% Function to analyze the best fit from the matrix of design evaluations, priority, and goals
function analyze_best_fit(P, Z, priority, weighted_priority, goals, design_constraints)
    % declare the indices
    high_priority_indices = [];
    moderate_priority_indices = [];
    low_priority_indices = [];

    % get the hard priority index
    hard_priority_index = find(priority == max(priority));

    if max(priority) - 1 >= 0
        % get the high priority index
        high_priority_indices = find(priority == max(priority) - 1);
    end

    if max(priority) - 2 >= 0
        % get the moderate priority index
        moderate_priority_indices = find(priority == max(priority) - 2);
    end

    if max(priority) - 3 >= 0
        % get the low priority index
        low_priority_indices = find(priority == max(priority) - 3);
    end

    % Get the indices that satisfy the hard priority
    hard_indices = find(Z(:, hard_priority_index) < goals(hard_priority_index));
    high_indices = [];
    moderate_indices = [];
    low_indices = [];

    % get the high priority index
    for i = 1:length(high_priority_indices)
        high_indices = [high_indices; find(Z(:, high_priority_indices(i)) < goals(high_priority_indices(i)))];
    end

    % get the moderate priority index
    for i = 1:length(moderate_priority_indices)
        moderate_indices = [moderate_indices; find(Z(:, moderate_priority_indices(i)) < goals(moderate_priority_indices(i)))];
    end

    % get the low priority index
    for i = 1:length(low_priority_indices)
        low_indices = [low_indices; find(Z(:, low_priority_indices(i)) < goals(low_priority_indices(i)))];
    end

    % get the intersection of the hard, high, moderate, and low indices
    hard_high_indices = intersect(hard_indices, high_indices);
    hard_high_moderate_indices = intersect(hard_high_indices, moderate_indices);
    hard_high_moderate_low_indices = intersect(hard_high_moderate_indices, low_indices);

    if isempty(hard_indices)
        fprintf('No indices satisfy the hard constraints, Bad design\n');
        indices = [];
    elseif isempty(hard_high_indices)
        fprintf('No indices satisfy the hard and high constraints\n');
        fprintf('The number of indices that satisfy the hard constraints is: %d\n', length(hard_indices));
        indices = hard_indices;
    elseif isempty(hard_high_moderate_indices)
        fprintf('No indices satisfy the hard, high, and moderate constraints\n');
        fprintf('The number of indices that satisfy the hard and high constraints is: %d\n', length(hard_high_indices));
        indices = hard_high_indices;
    elseif isempty(hard_high_moderate_low_indices)
        fprintf('No indices satisfy the hard, high, moderate, and low constraints\n');
        fprintf('The number of indices that satisfy the hard, high, and moderate constraints is: %d\n', length(hard_high_moderate_indices));
        indices = hard_high_moderate_indices;
    else
        fprintf('The number of indices that satisfy the hard, high, moderate, and low constraints is: %d\n', length(hard_high_moderate_low_indices));
        indices = hard_high_moderate_low_indices;
    end

    % check the difference between the goals and the Z values to find the best fit
    for idx = 1:length(indices)
        diff = abs(Z(indices(idx), :) - goals) / goals;
        diff = diff .* weighted_priority;
        % remove the NaN values
        diff(isnan(diff)) = 0;
        best_fit(idx) = sum(diff);
    end
    % get the index of the minimum best fit
    min_best_fit_index = find(best_fit == min(best_fit));

    % print the row of Z that satisfies the minimum best fit
    best_solution = Z(indices(min_best_fit_index), :);

    % print the violation of the constraints
    for i = 1:length(best_solution)
        if best_solution(i) > goals(i)
            % print the label of the feature
            fprintf('The %s constraint is violated\n', design_constraints{i});
        end
    end

    % print success rate deviation from goal
    success_rate = (abs(goals) - abs(best_solution)) ./ goals * 100;
    % round the success rate to 1 decimal place
    success_rate = round(success_rate, 1);

    % print the success rate
    fprintf('The success rate deviation from the goal is: %s\n', mat2str(success_rate));

    % print the best solution
    fprintf('The best solution is: %s\n', mat2str(best_solution));

    % get the values of Kp and Ki that satisfy the minimum best fit
    fprintf('The values of Kp and Ki that satisfy the minimum best fit are: %s\n', mat2str(P(indices(min_best_fit_index), :)));
end


% Function to build the optimizing engine with preferability
function buildOptimizingEngine(enable_preference, P, iterations, goals, priority, weighted_priority, best_sampling_plan, design_constraints)
    reference_point = max(optimizeControlSystem(P));
    convergence = zeros(1, iterations, 'double');
    bounds = [0, 0; 1, 1];
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    for i = 1:iterations
        % Step 1: Initializing the population
        Z = optimizeControlSystem(P);

        % Step 2: Calculating fitness
        % Step 2.1: Non-dominated sorting with preferability (flipped ranking)
        if enable_preference
            ranking = rank_prf(Z, goals, priority);
        else
            ranking = rank_nds(Z);
        end
        % inverse the ranking
        ranking = max(ranking) - ranking;

        % Step 2.2: Crowding distance assignment
        % NSGA-II density estimator
        distances = crowding(Z, ranking);

        % Step 3: Performing selection-for-selection
        % Binary tournament selection with replacement.
        % Returns the indices of the selected individuals.
        selectThese = btwr(distances, length(distances));

        % Step 4: Performing variation

        % Step 4.1: Simulated binary crossover
        % Simulated Binary Crossover operator
        % for real number representations.

        % Z -> objectives
        % P - > decision variables

        parents  = P(selectThese, :);
        offspring = sbx(parents, bounds);

        % Step 4.2: Polynomial mutation
        % Polynomial mutation operator
        % for real number representations.
        %
        C = polymut(offspring, bounds);

        % Step 5: Performing selection-for-survival

        % Step 5.1: Combine the parent and offspring populations
        unifiedPop = [P; C];

        % Step 5.2: Reducing the population
        % NSGA II clustering procedure.
        % Selects the new parent population from the unified population
        % of previous parents and offspring.

        Z_unified = optimizeControlSystem(unifiedPop);
        if enable_preference
            new_indices = reducerNSGA_II(unifiedPop, rank_prf(Z_unified, goals, priority), crowding(Z_unified, rank_prf(Z_unified, goals, priority)));
        else
            new_indices = reducerNSGA_II(unifiedPop, rank_nds(Z_unified), crowding(Z_unified, rank_nds(Z_unified)));
        end

        % Step 5.3: Select the new population
        P = unifiedPop(new_indices, :);

        % Step 6: Check for convergence
        convergence(i) = log10(Hypervolume_MEX(Z, reference_point));
        % Hypervolume is a measure of the volume of the objective space, dominated by the Pareto front.

        % Step 7: Plot the new population and convergence

        % Plot the progress using drawnow
        subplot(1, 2, 1);
        plot(P(:,1), P(:,2), 'o');
        title(sprintf('Sampling update for iteration %d', i));
        xlabel('x_1');
        ylabel('x_2');
        drawnow;
    end

    % Plot the convergence using subplot and lines
    subplot(1, 2, 2);
    plot(convergence, 'LineWidth', 2, 'Color', 'r');
    title('Hypervolume Convergence plot');
    xlabel('Iteration');
    ylabel('Hypervolume (log10 scale)');
    xlim([0, iterations + 1]);
    ylim([55, 56]);

    % Get the design evaluations
    Z = optimizeControlSystem(P);

    % Implement knowledge discovery
    knowledge_discovery(Z, best_sampling_plan, design_constraints, goals);

    % Analyze best fit from the matrix of design evaluations
    analyze_best_fit(P, Z, priority, weighted_priority, goals, design_constraints);
end

% Function to optimize the sampling plan
function Z_optimized = optimizeControlSystem(P)
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

        phi_metric = mmphi(P * scale, 5, 1);
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
    pause(1);
end

% Function to implement knowledge discovery
function knowledge_discovery(Z, best_sampling_plan, design_constraints, goals)
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    stability = nan * ones(size(Z, 2), 1);
    stability(1:3) = goals(1:3);

    transient = nan * ones(size(Z, 2), 1);
    transient(4:8) = goals(4:8);

    steady_state = nan * ones(size(Z, 2), 1);
    steady_state(9) = goals(9);

    sustainability = nan * ones(size(Z, 2), 1);
    sustainability(10) = goals(10);

    % Append extra line to A
    Z = [Z; stability'; transient'; steady_state'; sustainability'];

    % create a vector of strings
    groupDataVector = cell(1, size(Z, 1) - 4);
    for i = 1:length(groupDataVector)
        groupDataVector{i} = 'All Design Evaluation';
    end

    % Grouping vector
    groupDataVector = [groupDataVector, 'Stability', 'Transient', 'Steady State', 'Sustainability'];

    p = parallelplot(Z, 'groupData', groupDataVector, 'Color', {'green','red', 'blue', 'magenta', 'black'}, 'LineWidth', 2);
    p.MarkerStyle = {'none','o', 'o', 'o', 'o'};
    p.MarkerSize(end) = 15;
    p.LineStyle = {'-', '--', '--', '--', '--'};
    p.LegendTitle = 'Performance Criteria';


    p.CoordinateTickLabels = design_constraints;
    p.YLabel = 'Performance metric value';
    p.Title = sprintf('Performance evaluations for %s sampling plan', best_sampling_plan);
end

% Function to find the lower dimensional representation of the data
function data_mining(Z)
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    % remove inf values
    Z(isinf(Z)) = 1e6;

    % normalize the data
    Z = normalize(Z);

    % apply pca
    [coeff, score, latent, ~, explained] = pca(Z);

    % plot the explained variance
    subplot(2, 2, 1);
    plot(explained, 'o-');
    title('Explained Variance');
    xlabel('Principal Component');
    ylabel('Explained Variance (%)');

    % plot the first three principal components

    % perform kmeans clustering on the first three principal components
    [idx, C] = kmeans(score(:,1:3), 3);

    % plot the first three principal components
    subplot(2, 2, 3);
    scatter3(score(:,1), score(:,2), score(:,3), 50, idx, 'o');
    % change the color of the centroids
    hold on;
    scatter3(C(:,1), C(:,2), C(:,3), 100, 'k', 'x');
    hold off;
    title('First Three Principal Components');
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    % label the legend
    legend('Clusters', 'Centroids');

    [idx, C] = kmeans(score(:,4:6), 3);

    % plot the next three principal components
    subplot(2, 2, 4);
    scatter3(score(:,4), score(:,5), score(:,6), 50, idx, 'o');
    % change the color of the centroids
    hold on;
    scatter3(C(:,1), C(:,2), C(:,3), 100, 'k', 'x');
    hold off;
    title('Next Three Principal Components');
    xlabel('PC4');
    ylabel('PC5');
    zlabel('PC6');
    % label the legend
    legend('Clusters', 'Centroids');
    pause(0.5);
end
