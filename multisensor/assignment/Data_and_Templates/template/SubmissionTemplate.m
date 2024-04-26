function [x_est,b_est,Ax_f_instance,Ay_f_instance,Az_f_instance,p_f_instance,q_f_instance,r_f_instance,AoA_f_instance] = SID21010000(c_k, d_k, t, dt)

%% Note: Rename the function in the format of SID + your student ID as on Blackboard (e.g. if your ID is 21010000, name the function as SID21010000 and submit it as SID21010000.m)
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


end