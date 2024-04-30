function [x_est, b_est, Ax_f_instance, Ay_f_instance, Az_f_instance, p_f_instance, q_f_instance, r_f_instance, AoA_f_instance] = SID230118120(c_k, d_k, t, dt)
    %%%%%%%%%%%%%%%%%%%%%%%
    % Input Variables
    %   c_k: input measurements, a N x 6 matrix with rows the recordings of different time instances and the columns being [A_x_m, A_y_m, A_z_m, p_m, q_m, r_m]
    %   d_k: output measurements, a N x 12 matrix with rows the recordings of different time instances and the columns being [x_GPS_m, y_GPS_m, z_GPS_m, u_GPS_m, v_GPS_m, w_GPS_m, phi_GPS_m, theta_GPS_m, psi_GPS_m, V_TAS_m, alpha_m, beta_m]
    %   t: time vector
    %   dt: uniform time step size
    %%%%%%%%%%%%%%%%%%%%%%%
    % Output Variables
    %   x_est: estimated state trajectory, a N x 12 matrix with rows the recordings of different time instances and the columns being [x_E, y_E, z_E, u, v, w, \phi, \theta, \psi, V_{wxE}, V_{wyE}, V_{wzE}]
    %   b_est: estimated state trajectory, a N x 6 matrix with rows the recordings of different time instances and the columns being [b_A_x, b_A_y, b_A_z, b_p, b_q, b_r]
    %   Ax_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_x. In case of no faults detected, please return a value of 0,
    %   Ay_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_y. In case of no faults detected, please return a value of 0,
    %   Az_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of A_z. In case of no faults detected, please return a value of 0,
    %   p_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of p. In case of no faults detected, please return a value of 0,
    %   q_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of q. In case of no faults detected, please return a value of 0,
    %   r_f_instance: a vector with the values of corresponding time instances t for the detection of faults in measurements of r. In case of no faults detected, please return a value of 0,
    %   AoA_f_instance: a number equal to the corresponding time instances t for the first detection of faults in measurements of angle of attack sensor. In case of no faults detected, please return a value of 0.
    %%%%%%%%%%%%%%%%%%%%%%%

    %% Your code here

    % Define the data
    z_k = d_k; % output
    u_k = c_k; % input

    % Define the standard deviations
    sigma_A_x = 0.01;
    sigma_A_y = 0.01;
    sigma_A_z = 0.01;

    sigma_p = deg2rad(0.01);
    sigma_q = deg2rad(0.01);
    sigma_r = deg2rad(0.01);

    sigma_w_b_A_x = 1;
    sigma_w_b_A_y = 1;
    sigma_w_b_A_z = 1;

    sigma_w_b_p = deg2rad(1);
    sigma_w_b_q = deg2rad(1);
    sigma_w_b_r = deg2rad(1);

    sigma_x_E = 5;
    sigma_y_E = 5;
    sigma_z_E = 10;

    sigma_u = 0.1;
    sigma_v = 0.1;
    sigma_w = 0.1;
    sigma_V_tas = 0.1;

    sigma_phi = deg2rad(0.1);
    sigma_theta = deg2rad(0.1);
    sigma_psi = deg2rad(0.1);
    sigma_alpha = deg2rad(0.1);
    sigma_beta = deg2rad(0.1);

    % Define the initial conditions
    x_E = z_k(1,1); % x_GPS
    y_E = z_k(1,2); % y_GPS
    z_E = z_k(1,3); % z_GPS

    u_estimate = 85; % V_tas estimate
    v = 0;
    w = 0;

    b_A_x_estimate = 0.01;
    b_A_y_estimate = 0.01;
    b_A_z_estimate = 0.01;

    b_p_estimate = deg2rad(0.1);
    b_q_estimate = deg2rad(0.1);
    b_r_estimate = deg2rad(0.1);

    phi = z_k(1,7); % phi_GPS
    theta = z_k(1,8); % theta_GPS
    psi = z_k(1,9); % psi_GPS

    V_wxE = 0;
    V_wyE = 0;
    V_wzE = 0;

    % Define state names
    state_names = {'x_{E}', 'y_{E}', 'z_{E}', 'u', 'v', 'w', '\phi', '\theta', '\psi', 'b_{A_{x}}', 'b_{A_{y}}', 'b_{A_{z}}', 'b_{p}', 'b_{q}', 'b_{r}', 'V_{wxE}', 'V_{wyE}', 'V_{wzE}'};
    output_state_names = {'x_{GPS}', 'y_{GPS}', 'z_{GPS}', 'u_{GPS}', 'v_{GPS}', 'w_{GPS}', '\phi_{GPS}', '\theta_{GPS}', '\psi_{GPS}', 'V_{tas}', '\alpha', '\beta'};
    units = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad', 'm/s^2', 'm/s^2', 'm/s^2', 'rad/s', 'rad/s', 'rad/s', 'm/s', 'm/s', 'm/s'};
    output_units = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad', 'm/s', 'rad', 'rad'};

    stdw = [sigma_A_x sigma_A_y sigma_A_z sigma_p sigma_q sigma_r sigma_w_b_A_x sigma_w_b_A_y sigma_w_b_A_z sigma_w_b_p sigma_w_b_q sigma_w_b_r]; % standard deviation of process noise
    stdv = [sigma_x_E sigma_y_E sigma_z_E sigma_u sigma_v sigma_w sigma_phi sigma_theta sigma_psi sigma_V_tas sigma_alpha sigma_beta];      % standard deviation of measurement noise
    Ex_0 = [x_E y_E z_E u_estimate v w phi theta psi b_A_x_estimate b_A_y_estimate b_A_z_estimate b_p_estimate b_q_estimate b_r_estimate V_wxE V_wyE V_wzE]; % initial state estimate
    stdx_0 = [0.5 0.5 0.5 90 90 90 0.5 0.5 0.5 5 5 5 5 5 5 90 90 90];  % standard deviation of x_0

    % Run the Extended Kalman Filter
    [x_est, b_est, fault_instances_list] = ...
        runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0, state_names, output_state_names, units, output_units);

    % Extract the fault instances
    Ax_f_instance = fault_instances_list{1};
    Ay_f_instance = fault_instances_list{2};
    Az_f_instance = fault_instances_list{3};
    p_f_instance = fault_instances_list{4};
    q_f_instance = fault_instances_list{5};
    r_f_instance = fault_instances_list{6};
    AoA_f_instance = 0;
end


% Function to run the Extended Kalman Filter
function [x_est, b_est, fault_instances_list] = ...
        runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0, state_names, output_state_names, units, output_units)
    Ts = dt;     % time step (already provided by the data)
    N = length(t); % total number of steps

    xhat_km1_km1 = Ex_0; % x(0|0) = E{x_0}
    P_km1_km1 = diag(stdx_0.^2);  % P(0|0) = P(0)
    Q=diag(stdw.^2);
    R=diag(stdv.^2);

    n = length(xhat_km1_km1); % n: state dimension
    m = size(u_k, 2);     % m: observation dimension
    p = size(z_k, 2);     % m: observation dimension
    u_km1 = [zeros(1,m); u_k]; % shifted to have the right indices

    % Preallocate storage
    stdx_cor  = zeros(N, n);  % \sigma(k-1|k-1), standard deviation of state estimation error (hint: diagonal elements of P(k-1|k-1))
    x_cor     = zeros(N, n);  % \hat{x}(k-1|k-1), previous estimation
    K       = cell(N, 1);   % K(k) Kalman Gain
    innov     = zeros(N, p);  % y(k)-y(k|k-1), innovation, with y(k|k-1)=h(\hat{x}(k|k-1),u(k|k-1),k);

    for k=1:N
        % Step 1: Prediction
        [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
        xhat_k_km1 = x_nonlin(end,:); % x(k|k-1) (prediction)

        % Step 2: Covariance matrix of state prediction error / Minimum
        % prediction MSE
        [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts); % Phi(k,k-1), Gamma(k,k-1)
        P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1'; % P(k|k-1) (prediction)

        % Step 3: Kalman Gain
        H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
        Ve = (H_k * P_k_km1 * H_k' + R); % Pz(k|k-1) (prediction)
        K_k = P_k_km1 * H_k' / Ve; % K(k) (gain)

        % Step 4: Measurement Update (Correction)
        z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]); % z(k|k-1) (prediction of output)
        xhat_k_k = xhat_k_km1 + (z_k(k,:) - z_k_km1)*K_k'; % x(k|k) (correction)

        % Cyber attack on alpha, fix the measurement update of the 11th state of innovation
        % Innovation data calculation
        innov(k,:) = z_k(k,:) - z_k_km1; % y(k)-y(k|k-1) (innovation)

        % standardised the innovation
        innov(k, :) = innov(k, :) ./ sqrt(diag(Ve))';

        % Check if the innovation is within 3 standard deviations for the 11th state
        if abs(innov(k, 11)) > 3
            % If the innovation is outside 3 standard deviations, then the measurement is not used
            z_k(k, 11) = z_k_km1(11);
            innov(k, 11) = 0;
        end

        % Step 5: Correction for Covariance matrix of state Estimate error /
        % Minimum MSE
        I_KH = eye(n) - K_k * H_k;
        P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k'; % P(k|k) (correction)

        % Save data: State estimate and std dev
        stdx_cor(k,:) = sqrt(diag(P_km1_km1)); % \sigma(k-1|k-1) Standard deviation of state estimation error
        x_cor(k,:) = xhat_km1_km1; % \hat{x}(k-1|k-1), estimated state
        K{k,1} = K_k; % K(k) (gain)
        innov(k,:)= z_k(k,:) - z_k_km1;

        % Recursive step
        xhat_km1_km1 = xhat_k_k;
        P_km1_km1 = P_k_k;
    end
    % Return the estimated states
    x_est = [x_cor(:,1:9) x_cor(:,16:18)];
    b_est = x_cor(:,10:15);

    % Plot the estimated states
    plot_estimated_states(x_cor, t, N, innov, state_names, output_state_names, units, output_units);

    % Implement CUSUM test for fault detection
    x_faulty = b_est;
    theta0_list = mean(x_faulty);
    sigma0_list = std(x_faulty);
    thresholds_list = [6 0 1 0.2 0.1 0.1];  % empirical values
    cum_threshold_list = thresholds_list + theta0_list;
    data_labels = {'A_{x}', 'A_{y}', 'A_{z}', 'p', 'q', 'r'};
    leak = [1 1 2 1 1 1];

    % Implement CUSUM algorithm for states
    state_instances_list = implement_cusum(x_faulty, leak, cum_threshold_list, theta0_list, sigma0_list, data_labels);

    % Implement CUSUM algorithm for angle of attack
    x_faulty = innov(:,11);
    theta0 = mean(x_faulty);
    sigma0 = std(x_faulty);
    thresholds_list = [0.02];  % empirical values
    cum_threshold_list = thresholds_list + theta0;
    data_labels = {'\alpha'};
    leak = 1;

    alpha_instances = implement_cusum(x_faulty, leak, cum_threshold_list, theta0, sigma0, data_labels);

    % Save the fault instances
    fault_instances_list = [state_instances_list alpha_instances];
end

% Function to implement CUSUM algorithm
function fault_instances_list = implement_cusum(x_faulty, leak, cum_threshold_list, theta0_list, sigma0_list, data_labels)
    % Initialize variables
    skip_faulty_indices = 500; % Skip initial values so that EKF can converge
    faulty_size = size(x_faulty, 1);
    faulty_columns = size(x_faulty, 2);
    fault_instances = cell(1, faulty_columns);

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

        % Add the skipped indices
        if ~isempty(k_alarm_pos)
            k_alarm_pos = k_alarm_pos + skip_faulty_indices;
        end
        if ~isempty(k_alarm_neg)
            k_alarm_neg = k_alarm_neg + skip_faulty_indices;
        end

        % Save the fault instances
        if ~isempty(k_alarm_pos) && ~isempty(k_alarm_neg)
            % Find the union of the two alarms
            fault_instances{idx} = unique([k_alarm_pos; k_alarm_neg]);
        elseif ~isempty(k_alarm_pos)
            fault_instances{idx} = k_alarm_pos;
        elseif ~isempty(k_alarm_neg)
            fault_instances{idx} = k_alarm_neg;
        else
            fault_instances{idx} = 0;
        end
        fault_instances{idx} = fault_instances{idx} / 100; % Convert to seconds

        figure
        p1=plot(g_pos);
        hold on
        p2=plot(g_neg);

        % Save the fault instances
        if ~isempty(k_alarm_pos) && ~isempty(k_alarm_neg)

            for i=1:length(k_alarm_pos)
                p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-20 20],'b--');
            end
            for i=1:length(k_alarm_neg)
                p4=plot([k_alarm_neg(i) k_alarm_neg(i)],[-20 20],'r-.');
            end
            legend([p1,p2,p3,p4],'Positive Test', 'Negative Test','Alarms for positive test','Alarms for negative test')

        elseif ~isempty(k_alarm_pos)

            for i=1:length(k_alarm_pos)
                p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-20 20],'b--');
            end
            legend([p1,p2,p3],'Positive Test', 'Negative Test','Alarms for positive test')

        elseif ~isempty(k_alarm_neg)

            for i=1:length(k_alarm_neg)
                p3=plot([k_alarm_neg(i) k_alarm_neg(i)],[-20 20],'r-.');
            end
            legend([p1,p2,p3],'Positive Test', 'Negative Test','Alarms for negative test')

        else
            legend([p1,p2],'Positive Test', 'Negative Test')
        end

        yline(threshold_neg)
        yline(threshold_pos)
        ylim([threshold_neg-0.025 threshold_pos+0.025])
        xlabel('Step')
        ylabel('g_t')
        title(['CUSUM Algorithm for ', data_labels{idx}]);
    end

    fault_instances_list = fault_instances;
end


function plot_estimated_states(x_cor, t, N, innov, state_names, output_state_names, units, output_units)
    % Plot the estimated states
    x_label = 'Time [s]';
    y_label = 'Estimation in ';
    font_size = 12;
    n_rows = ceil(length(state_names) / 3);

    figure
    for i = 1:length(state_names)
        subplot(n_rows, 3, i)
        plot(t(1:N), x_cor(1:N,i), '--b', 'LineWidth', 2)
        grid on
        xlabel(x_label, 'FontSize', font_size)
        ylabel(y_label + " (" + units(i) + ")", 'FontSize', font_size)
        title(state_names(i), 'FontSize', font_size)
    end
    % set figure name
    set(gcf, 'Name', 'Estimated States through Extended Kalman Filter with Biases and Faults')

    % Plot the innovation
    figure
    for i = 1:size(innov, 2)
        subplot(n_rows, 3, i)
        plot(t(1:N), innov(1:N,i), '--r', 'LineWidth', 2)
        grid on
        xlabel(x_label, 'FontSize', font_size)
        ylabel('Innovation' + " (" + output_units(i) + ")", 'FontSize', font_size)
        title(output_state_names(i), 'FontSize', font_size)
    end
    % set figure name
    set(gcf, 'Name', 'Innovation through Extended Kalman Filter with Faults')


    % Plot bias
    figure
    for i = 10:15
        subplot(2, 3, i-9)
        plot(t(1:N), x_cor(1:N,i), '--b', 'LineWidth', 2)
        grid on
        xlabel(x_label, 'FontSize', font_size)
        ylabel(y_label + " (" + units(i) + ")", 'FontSize', font_size)
        title(state_names(i), 'FontSize', font_size)
    end
end

% Function to calculate the state vector derivative
function x_dot_vector = funcf(x_vector, c_vector, t)
    % Extract values from the state vector
    x_E = x_vector(1);
    y_E = x_vector(2);
    z_E = x_vector(3);

    u = x_vector(4);
    v = x_vector(5);
    w = x_vector(6);

    phi = x_vector(7);
    theta = x_vector(8);
    psi = x_vector(9);

    b_A_x = x_vector(10);
    b_A_y = x_vector(11);
    b_A_z = x_vector(12);

    b_p = x_vector(13);
    b_q = x_vector(14);
    b_r = x_vector(15);

    V_wxE = x_vector(16);
    V_wyE = x_vector(17);
    V_wzE = x_vector(18);

    % Extract values from the input vector
    A_x = c_vector(1);
    A_y = c_vector(2);
    A_z = c_vector(3);

    p = c_vector(4);
    q = c_vector(5);
    r = c_vector(6);

    % Gravity
    g = 9.81;

    % Kinematic Model with Biases and Faults (excluding noises)
    x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+V_wxE;
    y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+V_wyE;
    z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+V_wzE;
    u_dot =(A_x-b_A_x)-g*sin(theta)+(r-b_r)*v-(q-b_q)*w;
    v_dot =(A_y-b_A_y)+g*cos(theta)*sin(phi)+(p-b_p)*w-(r-b_r)*u;
    w_dot =(A_z-b_A_z)+g*cos(theta)*cos(phi)+(q-b_q)*u-(p-b_p)*v;
    phi_dot =(p-b_p)+(q-b_q)*sin(phi)*tan(theta)+(r-b_r)*cos(phi)*tan(theta);
    theta_dot =(q-b_q)*cos(phi)-(r-b_r)*sin(phi);
    psi_dot =(q-b_q)*sin(phi)/cos(theta)+(r-b_r)*cos(phi)/cos(theta);
    b_A_x_dot =0;
    b_A_y_dot =0;
    b_A_z_dot =0;
    b_p_dot =0;
    b_q_dot =0;
    b_r_dot =0;
    V_wxE_dot =0;
    V_wyE_dot =0;
    V_wzE_dot =0;

    x_dot_vector = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot b_A_x_dot b_A_y_dot b_A_z_dot b_p_dot b_q_dot b_r_dot V_wxE_dot V_wyE_dot V_wzE_dot]';
end

% Function to calculate the observation vector
function d_vector = funch(x_vector, c_vector, t)
    % Extract values from the state vector
    x_E = x_vector(1);
    y_E = x_vector(2);
    z_E = x_vector(3);

    u = x_vector(4);
    v = x_vector(5);
    w = x_vector(6);

    phi = x_vector(7);
    theta = x_vector(8);
    psi = x_vector(9);

    b_A_x = x_vector(10);
    b_A_y = x_vector(11);
    b_A_z = x_vector(12);

    b_p = x_vector(13);
    b_q = x_vector(14);
    b_r = x_vector(15);

    V_wxE = x_vector(16);
    V_wyE = x_vector(17);
    V_wzE = x_vector(18);

    % Observation Model [Output]
    x_GPS = x_E;
    y_GPS = y_E;
    z_GPS = z_E;
    u_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+V_wxE;
    v_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+V_wyE;
    w_GPS = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+V_wzE;
    phi_GPS = phi;
    theta_GPS = theta;
    psi_GPS = psi;
    V_tas = sqrt(u^2+v^2+w^2);
    alpha = atan(w/u);
    beta = atan(v/sqrt(u^2+w^2));

    d_vector = [x_GPS y_GPS z_GPS u_GPS v_GPS w_GPS phi_GPS theta_GPS psi_GPS V_tas alpha beta];
end

% Function to linearize the discrete-time dynamics
function [Phi, Gamma] = funcLinDisDyn(x_vector, c_m_vector, Ts)

    % Extract values from the state vector
    x_E = x_vector(1);
    y_E = x_vector(2);
    z_E = x_vector(3);

    u = x_vector(4);
    v = x_vector(5);
    w = x_vector(6);

    phi = x_vector(7);
    theta = x_vector(8);
    psi = x_vector(9);

    b_A_x = x_vector(10);
    b_A_y = x_vector(11);
    b_A_z = x_vector(12);

    b_p = x_vector(13);
    b_q = x_vector(14);
    b_r = x_vector(15);

    V_wxE = x_vector(16);
    V_wyE = x_vector(17);
    V_wzE = x_vector(18);

    % Extract values from the input measurement vector
    A_x_m = c_m_vector(1);
    A_y_m = c_m_vector(2);
    A_z_m = c_m_vector(3);
    p_m = c_m_vector(4);
    q_m = c_m_vector(5);
    r_m = c_m_vector(6);

    % Numerical evaluation of continuous - time dynamics
    % input measurement noises

    A = [0, 0, 0, cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)),                                   -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)),  0,  0,  0,  0,                    0,                    0, 1, 0, 0;
         0, 0, 0, cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)),                                   -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)),  0,  0,  0,  0,                    0,                    0, 0, 1, 0;
         0, 0, 0,         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),                                             - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 1;
         0, 0, 0,                   0,                                        r_m - b_r,                                        b_q - q_m,                                                                                  0,                                                                             -(981*cos(theta))/100,                                                                                                     0, -1,  0,  0,  0,                    w,                   -v, 0, 0, 0;
         0, 0, 0,           b_r - r_m,                                                0,                                        p_m - b_p,                                                      (981*cos(phi)*cos(theta))/100,                                                                    -(981*sin(phi)*sin(theta))/100,                                                                                                     0,  0, -1,  0, -w,                    0,                    u, 0, 0, 0;
         0, 0, 0,           q_m - b_q,                                        b_p - p_m,                                                0,                                                     -(981*cos(theta)*sin(phi))/100,                                                                    -(981*cos(phi)*sin(theta))/100,                                                                                                     0,  0,  0, -1,  v,                   -u,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                  sin(phi)*tan(theta)*(b_r - r_m) - cos(phi)*tan(theta)*(b_q - q_m),               - cos(phi)*(b_r - r_m)*(tan(theta)^2 + 1) - sin(phi)*(b_q - q_m)*(tan(theta)^2 + 1),                                                                                                     0,  0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta), 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                        cos(phi)*(b_r - r_m) + sin(phi)*(b_q - q_m),                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,            -cos(phi),             sin(phi), 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,              (sin(phi)*(b_r - r_m))/cos(theta) - (cos(phi)*(b_q - q_m))/cos(theta), - (cos(phi)*sin(theta)*(b_r - r_m))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_q - q_m))/cos(theta)^2,                                                                                                     0,  0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta), 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0];


    G = [0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         -1,  0,  0,  0,                    w,                   -v,    1,    0,    0,  0,                  -w,                   v;
         0, -1,  0, -w,                    0,                    u,    0,    1,    0,  w,                   0,                  -u;
         0,  0, -1,  v,                   -u,                    0,    0,    0,    1, -v,                   u,                   0;
         0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),    0,    0,    0,  1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0,  0,  0,  0,            -cos(phi),             sin(phi),    0,    0,    0,  0,            cos(phi),           -sin(phi);
         0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),    0,    0,    0,  0, sin(phi)/cos(theta), cos(phi)/cos(theta);
         0,  0,  0,  0,                    0,                    0, 1/Ts,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0, 1/Ts,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0, 1/Ts,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,  0,                   0,                   0];


    G = [0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
         -1,  0,  0,  0,                    w,                   -v,    1,    0,    0,    0,                  -w,                   v;
         0, -1,  0, -w,                    0,                    u,    0,    1,    0,    w,                   0,                  -u;
         0,  0, -1,  v,                   -u,                    0,    0,    0,    1,   -v,                   u,                   0;
         0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),    0,    0,    0,    1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0,  0,  0,  0,            -cos(phi),             sin(phi),    0,    0,    0,    0,            cos(phi),           -sin(phi);
         0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),    0,    0,    0,    0, sin(phi)/cos(theta), cos(phi)/cos(theta);
         0,  0,  0,  0,                    0,                    0, 1/Ts,    0,    0,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0, 1/Ts,    0,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0, 1/Ts,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0, 1/Ts,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                1/Ts,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                1/Ts;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0;
         0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,                   0,                   0];


    % Discretisation of dynamics
    [Phi, Gamma] = c2d(A, G, Ts);

end

% Function to linearize the observation model
function H = funcLinDisObs(x_vector, c_m_vector, t)
    % Extract values from state vector
    x_E = x_vector(1);
    y_E = x_vector(2);
    z_E = x_vector(3);

    u = x_vector(4);
    v = x_vector(5);
    w = x_vector(6);

    phi = x_vector(7);
    theta = x_vector(8);
    psi = x_vector(9);

    b_A_x = x_vector(10);
    b_A_y = x_vector(11);
    b_A_z = x_vector(12);

    b_p = x_vector(13);
    b_q = x_vector(14);
    b_r = x_vector(15);

    V_wxE = x_vector(16);
    V_wyE = x_vector(17);
    V_wzE = x_vector(18);

    % Extract values from input measurement vector
    A_x_m = c_m_vector(1);
    A_y_m = c_m_vector(2);
    A_z_m = c_m_vector(3);
    p_m = c_m_vector(4);
    q_m = c_m_vector(5);
    r_m = c_m_vector(6);

    % Numerical evaluation of continuous - time dynamics
    H = [1, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 1, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 1,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0,                              cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 0, 0, 0, 0, 0, 0, 1, 0, 0;
         0, 0, 0,                              cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)), -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 0, 0, 0, 0, 0, 0, 1, 0;
         0, 0, 0,                                      -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
         0, 0, 0,                                                0,                                                0,                                                0,                                                                                  1,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               1,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                        v/(u^2 + v^2 + w^2)^(1/2),                        w/(u^2 + v^2 + w^2)^(1/2),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                0,                              1/(u*(w^2/u^2 + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),      1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)), -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
end
