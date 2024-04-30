% Copyright (c) 2024 Leander Stephen D'Souza

% Program to implement CUSUM algorithm for fault detection

% Read data files
data_file_nominal = load('saved_task2_2.mat');
data_file_faulty = load('saved_task3_1.mat');

% Define data labels
data_labels = {'A_{x}', 'A_{y}', 'A_{z}', 'p', 'q', 'r'};

% Extract required data
x_nominal = data_file_nominal.x_cor(:,10:15);
x_faulty = data_file_faulty.x_cor(:,10:15);

% Precomputed values for nominal data
theta0_list = [0.5000 0.3009 -0.0950 0.0052 -0.0070 0.0105]; % mean of the noise
sigma0_list = [0.0010 0.0021 0.0000 0.0000 0.0000 0.0000]; % variance of the noise
thresholds_list = [6 0 1 0.2 0.1 0.1];  % empirical values
cum_threshold_list = thresholds_list + theta0_list;


% Define leakage values
leak = [1 1 2 1 1 1];


% Implement CUSUM algorithm
implement_cusum(x_faulty, leak, cum_threshold_list, theta0_list, sigma0_list, data_labels);


% Function to implement CUSUM algorithm
function implement_cusum(x_faulty, leak, cum_threshold_list, theta0_list, sigma0_list, data_labels)
    % Initialize variables
    skip_faulty_indices = 500; % Skip initial values so that EKF can converge
    faulty_size = size(x_faulty, 1);
    faulty_columns = size(x_faulty, 2);

    for idx=1:faulty_columns
        straingauge = x_faulty(skip_faulty_indices:faulty_size, idx);
        theta0 = theta0_list(idx);
        sigma0 = sigma0_list(idx);
        cum_threshold = cum_threshold_list(idx);

        threshold_pos = cum_threshold;
        threshold_neg = -cum_threshold;

        % Two-Sided CUSUM Test
        g_pos = 0 * straingauge;
        g_neg = 0 * straingauge;
        k_alarm_pos = [];
        k_alarm_neg = [];

        s = (straingauge - theta0) / sigma0;

        for k = 1:size(straingauge, 1) - 1
            g_pos(k+1) = g_pos(k) + s(k) - leak(idx);
            g_neg(k+1) = g_neg(k) + s(k) + leak(idx);

            % Positive test
            if g_pos(k+1) < 0
                g_pos(k+1) = 0;
            end
            if g_pos(k+1) > threshold_pos
                k_alarm_pos = [k_alarm_pos; k+1];
                g_pos(k+1) = 0; % reset
            end

            % Negative test
            if g_neg(k+1) > 0
                g_neg(k+1) = 0;
            end
            if g_neg(k+1) < threshold_neg
                k_alarm_neg = [k_alarm_neg; k+1];
                g_neg(k+1) = 0; % reset
            end
        end

        % figure
        % p1=plot(g_pos);
        % hold on
        % p2=plot(g_neg);
        % for i=1:length(k_alarm_pos)
        %     p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-20 20],'b--');
        % end
        % for i=1:length(k_alarm_neg)
        %     p4=plot([k_alarm_neg(i) k_alarm_neg(i)],[-20 20],'r-.');
        % end
        % legend([p1,p2,p3,p4],'Positive Test', 'Negative Test','Alarms for positive test','Alarms for negative test')
        % yline(threshold_neg)
        % yline(threshold_pos)
        % ylim([threshold_neg-0.025 threshold_pos+0.025])
        % xlabel('Step')
        % ylabel('g_t')
        % title(['CUSUM Algorithm for ', data_labels{idx}]);
    end
end
