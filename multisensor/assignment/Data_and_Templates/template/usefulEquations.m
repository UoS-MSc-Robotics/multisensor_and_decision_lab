    %Bias Free Kinematic Model (excluding noises)
    x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
    y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
    z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    u_dot =A_x-g*sin(theta)+r*v-q*w;
    v_dot =A_y+g*cos(theta)*sin(phi)+p*w-r*u;
    w_dot =A_z+g*cos(theta)*cos(phi)+q*u-p*v;
    phi_dot =p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
    theta_dot =q*cos(phi)-r*sin(phi);
    psi_dot =q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta);
    v_wxE_dot =0;
    v_wyE_dot =0;
    v_wzE_dot =0;

    %Kinematic Model with Biases and Faults (excluding noises)
    x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
    y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
    z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
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
    v_wxE_dot =0;
    v_wyE_dot =0;
    v_wzE_dot =0;

    % Observation Model (excluding noises) (output)
    x_GPS=x;
    y_GPS=y;
    z_GPS=z;
    u_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
    v_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
    w_GPS=-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    phi_GPS=phi;
    theta_GPS=theta;
    psi_GPS=psi;
    V=sqrt(u^2+v^2+w^2);
    alpha=atan(w/u);
    beta=atan(v/sqrt(u^2+w^2));


