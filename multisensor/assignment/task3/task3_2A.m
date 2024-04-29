% Copyright (c) 2024 Leander Stephen D'Souza

% Program to linearize the system dynamics using symbolic toolbox with bias and faults

% Define symbolic variables

% IMU Acceleration
syms A_x A_y A_z
% IMU Angular Rates
syms p q r

% IMU Process noise
syms w_A_x w_A_y w_A_z w_p w_q w_r
% IMU Measurement noise
syms A_x_m A_y_m A_z_m p_m q_m r_m
% IMU bias for Measurement noise
syms b_A_x b_A_y b_A_z b_p b_q b_r
% IMU fault for Measurement noise
syms w_b_A_x w_b_A_y w_b_A_z w_b_p w_b_q w_b_r

% GPS Position
syms x_E y_E z_E

% Ground Earth speed
syms u v w
% Ground Earth angles
syms phi theta psi

% True air speed
syms V_tas
% Wind speed in Earth frame
syms V_wxE V_wyE V_wzE

% Angle of attack
syms alpha
% Sideslip angle
syms beta
% Sampling time
syms Ts


% Original continuous-time nonlinear system equations
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
b_A_x_dot = 1/Ts*(w_b_A_x);
b_A_y_dot = 1/Ts*(w_b_A_y);
b_A_z_dot = 1/Ts*(w_b_A_z);
b_p_dot = 1/Ts*(w_b_p);
b_q_dot = 1/Ts*(w_b_q);
b_r_dot = 1/Ts*(w_b_r);
V_wxE_dot =0;
V_wyE_dot =0;
V_wzE_dot =0;


% Observation Model (excluding noises) [Output]
x_GPS=x_E;
y_GPS=y_E;
z_GPS=z_E;
u_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+V_wxE;
v_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+V_wyE;
w_GPS=-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+V_wzE;
phi_GPS=phi;
theta_GPS=theta;
psi_GPS=psi;
V_tas=sqrt(u^2+v^2+w^2);
alpha=atan(w/u);
beta=atan(v/sqrt(u^2+w^2));

% Linearisation of the system

% State vector
x_vector = [x_E y_E z_E u v w phi theta psi b_A_x b_A_y b_A_z b_p b_q b_r V_wxE V_wyE V_wzE];

% Input vector
c_vector = [A_x A_y A_z p q r];

% Input bias vector
b_vector = [b_A_x b_A_y b_A_z b_p b_q b_r];

% Measurement vector
w_vector = [w_A_x w_A_y w_A_z w_p w_q w_r];
w_b_vector = [w_b_A_x w_b_A_y w_b_A_z w_b_p w_b_q w_b_r];

% Input measurement vector
c_m_vector = [A_x_m A_y_m A_z_m p_m q_m r_m];

% Output vector
d_vector = [x_GPS y_GPS z_GPS u_GPS v_GPS w_GPS phi_GPS theta_GPS psi_GPS V_tas alpha beta];

% Derivative of the state vector
x_dot_vector = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot b_A_x_dot b_A_y_dot b_A_z_dot b_p_dot b_q_dot b_r_dot V_wxE_dot V_wyE_dot V_wzE_dot];

% Substitute the input measurement equation in the system equations
x_dot_vector = subs(x_dot_vector, c_vector, c_m_vector - w_vector + w_b_vector);

% Calculate the Jacobian linearization
f_expression = subs(x_dot_vector, w_vector, zeros(1, length(w_vector)));
A_expression = jacobian(f_expression, x_vector);
G_expression = jacobian(x_dot_vector, [w_vector w_b_vector]);
H_expression = jacobian(d_vector, x_vector);

% Print the expressions
% disp(A_expression);
disp(G_expression);
% disp(H_expression);
