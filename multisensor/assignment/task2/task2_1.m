% Copyright (c) 2024 Leander Stephen D'Souza

% Program to implement Extended Kalman Filter for Aircraft Climb with Sensor Bias and Faults

% Read data files
datasets = [load('dataTask1.mat') load('dataTask2.mat')];
dataset_names = {'Unbiased Dataset', 'Biased Dataset'};

% Define intital values
x_cor_list = cell(length(datasets), 1);

for idx = 1:length(datasets)
    data_file = datasets(idx);

    % Define the data
    z_k = data_file.d_k; % output
    u_k = data_file.c_k; % input
    dt = data_file.dt; % time step
    t = data_file.t; % time vector

    % Define the standard deviations
    sigma_A_x = 0.01;
    sigma_A_y = 0.01;
    sigma_A_z = 0.01;

    sigma_p = deg2rad(0.01);
    sigma_q = deg2rad(0.01);
    sigma_r = deg2rad(0.01);

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

    u_estimate = 92
    v = 0;
    w = 0;

    phi = z_k(1,7); % phi_GPS
    theta = z_k(1,8); % theta_GPS
    psi = z_k(1,9); % psi_GPS

    V_wxE = 0;
    V_wyE = 0;
    V_wzE = 0;

    % Define state names
    state_names = {'x_{E}', 'y_{E}', 'z_{E}', 'u', 'v', 'w', '\phi', '\theta', '\psi', 'V_{wxE}', 'V_{wyE}', 'V_{wzE}'};
    units = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad', 'm/s', 'm/s', 'm/s'};

    stdw = [sigma_A_x sigma_A_y sigma_A_z sigma_p sigma_q sigma_r];     % standard deviation of system noise
    stdv = [sigma_x_E sigma_y_E sigma_z_E sigma_u sigma_v sigma_w sigma_phi sigma_theta sigma_psi sigma_V_tas sigma_alpha sigma_beta];      % standard deviation of measurement noise
    Ex_0 = [x_E y_E z_E u_estimate v w phi theta psi V_wxE V_wyE V_wzE]; % expected value of x_0
    stdx_0 = [0.5 0.5 0.5 90 90 90 0.5 0.5 0.5 90 90 90];  %standard deviation of x_0

    % Run the Extended Kalman Filter
    x_cor_list{idx} = runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0);
end

% Plot the estimated states
plot_estimated_states(x_cor_list, t, state_names, units, dataset_names);



% Function to run the Extended Kalman Filter
function x_cor = runEKF(u_k, z_k, t, dt, stdw, stdv, stdx_0, Ex_0)
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
end

function plot_estimated_states(x_cor_list, t, state_names, units, dataset_names)
    % Plot the estimated states
    x_label = 'Time [s]';
    y_label = 'Estimation in ';
    font_size = 12;

    figure
    for i = 1:length(state_names)
        subplot(4,3,i)
        hold on
        for idx = 1:length(x_cor_list)
            plot(t, x_cor_list{idx}(:,i), 'DisplayName', dataset_names{idx})
        end
        hold off
        grid on
        xlabel(x_label, 'FontSize', font_size)
        ylabel([y_label units{i}], 'FontSize', font_size)
        title([state_names{i}], 'FontSize', font_size)
        legend('Location', 'best')
    end
    % set figure name
    set(gcf, 'Name', 'Estimated States through Extended Kalman Filter with Sensor Bias and Faults')
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

    V_wxE = x_vector(10);
    V_wyE = x_vector(11);
    V_wzE = x_vector(12);

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
    u_dot =A_x-g*sin(theta)+r*v-q*w;
    v_dot =A_y+g*cos(theta)*sin(phi)+p*w-r*u;
    w_dot =A_z+g*cos(theta)*cos(phi)+q*u-p*v;
    phi_dot =p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
    theta_dot =q*cos(phi)-r*sin(phi);
    psi_dot =q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta);
    V_wxE_dot=0;
    V_wyE_dot=0;
    V_wzE_dot=0;

    x_dot_vector = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot V_wxE_dot V_wyE_dot V_wzE_dot]';
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

    V_wxE = x_vector(10);
    V_wyE = x_vector(11);
    V_wzE = x_vector(12);

    % Observation Model (excluding noises) [Output]
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

    V_wxE = x_vector(10);
    V_wyE = x_vector(11);
    V_wzE = x_vector(12);

    % Extract values from the input measurement vector
    A_x_m = c_m_vector(1);
    A_y_m = c_m_vector(2);
    A_z_m = c_m_vector(3);
    p_m = c_m_vector(4);
    q_m = c_m_vector(5);
    r_m = c_m_vector(6);

    % Numerical evaluation of continuous - time dynamics
    % input measurement noises

    A = [0, 0, 0, cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)),                 -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0;
         0, 0, 0, cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)),                 -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0;
         0, 0, 0,         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),                           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1;
         0, 0, 0,                   0,                                              r_m,                                             -q_m,                                                                                  0,                                                           -(981*cos(theta))/100,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                -r_m,                                                0,                                              p_m,                                                      (981*cos(phi)*cos(theta))/100,                                                  -(981*sin(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                 q_m,                                             -p_m,                                                0,                                                     -(981*cos(theta)*sin(phi))/100,                                                  -(981*cos(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                  q_m*cos(phi)*tan(theta) - r_m*sin(phi)*tan(theta),               r_m*cos(phi)*(tan(theta)^2 + 1) + q_m*sin(phi)*(tan(theta)^2 + 1),                                                                                                     0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                      - r_m*cos(phi) - q_m*sin(phi),                                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                              (q_m*cos(phi))/cos(theta) - (r_m*sin(phi))/cos(theta), (r_m*cos(phi)*sin(theta))/cos(theta)^2 + (q_m*sin(phi)*sin(theta))/cos(theta)^2,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                               0,                                                                                                     0, 0, 0, 0];

    G = [0,  0,  0,  0,                    0,                    0;
         0,  0,  0,  0,                    0,                    0;
         0,  0,  0,  0,                    0,                    0;
         -1,  0,  0,  0,                    w,                   -v;
         0, -1,  0, -w,                    0,                    u;
         0,  0, -1,  v,                   -u,                    0;
         0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta);
         0,  0,  0,  0,            -cos(phi),             sin(phi);
         0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta);
         0,  0,  0,  0,                    0,                    0;
         0,  0,  0,  0,                    0,                    0;
         0,  0,  0,  0,                    0,                    0];

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
    V_wxE = x_vector(10);
    V_wyE = x_vector(11);
    V_wzE = x_vector(12);

    % Extract values from input measurement vector
    A_x_m = c_m_vector(1);
    A_y_m = c_m_vector(2);
    A_z_m = c_m_vector(3);
    p_m = c_m_vector(4);
    q_m = c_m_vector(5);
    r_m = c_m_vector(6);

    % Numerical evaluation of continuous - time dynamics
    H = [1, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0;
         0, 1, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 1,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                              cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0;
         0, 0, 0,                              cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)), -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0;
         0, 0, 0,                                      -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1;
         0, 0, 0,                                                0,                                                0,                                                0,                                                                                  1,                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               1,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     1, 0, 0, 0;
         0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                        v/(u^2 + v^2 + w^2)^(1/2),                        w/(u^2 + v^2 + w^2)^(1/2),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                0,                              1/(u*(w^2/u^2 + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0;
         0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),      1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)), -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0];

end
