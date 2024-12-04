%% Fine TLC (Inbuilt)

clc; clear; close all;


t_ode = 0:1:1800;
tau_list = linspace(0.1, 120);

Jxf = zeros(size(tau_list));
Jxs = zeros(size(tau_list));
Juf = zeros(size(tau_list));
Jus = zeros(size(tau_list));

c = 0;
Tmax = 0.2;

for j = 1:1:length(tau_list)

tau = tau_list(j);

[T_final, yaw_d, yaw_found1, pitch_d, pitch_found1, roll_d, roll_found1] = fine_TLC(t_ode, tau, c);
[T_final_stan, yaw_d_stan, yaw_found1_stan, pitch_d_stan, pitch_found1_stan, roll_d_stan, roll_found1_stan] = standard_TLC(t_ode, tau, c);
Eul_fine = [roll_found1-roll_d; pitch_found1 - pitch_d; yaw_found1 - yaw_d];
Eul_stan = [roll_found1_stan-roll_d_stan; pitch_found1_stan - pitch_d_stan; yaw_found1_stan - yaw_d_stan];


T_mag = zeros(size(t_ode));
T_mag_stan = zeros(size(t_ode));
Eul_mag = zeros(size(t_ode));
Eul_mag_stan = zeros(size(t_ode));

for i=1:1:length(t_ode)
    T_mag(i) = T_final(:, i)'*T_final(:, i);
    T_mag_stan(i) = T_final_stan(:, i)'*T_final_stan(:, i);
    Eul_mag(i) = Eul_fine(:, i)'*Eul_fine(:, i);
    Eul_mag_stan(i) = Eul_stan(:, i)'*Eul_stan(:, i);
end

Ju_fine = trapz(t_ode, T_mag)/(2*Tmax*Tmax);
Ju_stan = trapz(t_ode, T_mag_stan)/(2*Tmax*Tmax);

Jx_fine = trapz(t_ode, Eul_mag)/(2*180*180);
Jx_stan = trapz(t_ode, Eul_mag_stan)/(2*180*180);

Jxf(j) = Jx_fine;
Jxs(j) = Jx_stan;
Juf(j) = Ju_fine;
Jus(j) = Ju_stan;

end

%% Plotting Block:

figure(4)
plot(tau_list, Jxf, 'LineWidth', 2, 'LineStyle','-', 'Color', 'b', 'DisplayName','Fine TLC')
hold on
plot(tau_list, Jxs, 'LineWidth', 2, 'LineStyle','-', 'Color', 'r', 'DisplayName','Traditional TLC')
legend("show")
ax = gca;
ax.FontSize = 16;
%ylim(ax, [0, 0.35])
grid on
title('Performance Index J_x Variation')
ylabel("Performance Index J_x -->")
xlabel("\tau (s) -->")

figure(5)
plot(tau_list, Juf, 'LineWidth', 2, 'LineStyle','-', 'Color', 'b', 'DisplayName','Fine TLC')
hold on
plot(tau_list, Jus, 'LineWidth', 2, 'LineStyle','-', 'Color', 'r', 'DisplayName','Traditional TLC')
legend("show")
ax = gca;
ax.FontSize = 16;
%ylim(ax, [0, 0.35])
grid on
title('Performance Index J_u Variation')
ylabel("Performance Index J_u -->")
xlabel("\tau (s) -->")


%%

function [T_final, yaw_d, yaw_found1, pitch_d, pitch_found1, roll_d, roll_found1] = fine_TLC(t_ode, tau, c)
%% Spacecraft model parameters (Wu et al., 2018):

t_list = t_ode;
dt = t_list(2)-t_list(1);

Js = (1+c)*[147, 6.5, 6; 6.5, 158, 5.5; 6, 5.5, 137];
Td = [1+2*sin(0.005*t_list); -1-5*cos(0.005*t_list); 2-4*cos(0.005*t_list)];

w0 = [1; -1; 1]*pi/180;
T0 = [0; 0; 0];

roll0 = 42*pi/180;
pitch0 = -19*pi/180;
yaw0 = -7*pi/180;

q0 = eul2quat([yaw0, pitch0, roll0])';
q0 = [q0(2); q0(3); q0(4); q0(1)];

Tmax = 0.2;

%% Desired trajectory:

if(t_list(end)>1100)

    t1100_ind = find(t_list==1100);
    
    t1 = t_list(1:t1100_ind-1);
    t2 = t_list(t1100_ind:end);
    
    roll_d = [-90*cos(0.025*t1), 0*t2]*pi/180;
    pitch_d = [-30*sin(0.005*t1), 0*t2]*pi/180;
    yaw_d = [75*ones(1, length(t1)), 0*t2]*pi/180;

else

    roll_d = -90*cos(0.025*t_list)*pi/180;
    pitch_d = -30*sin(0.005*t_list)*pi/180;
    yaw_d = 75*ones(1, length(t_list))*pi/180;

end

w_kin = 0.01;
w_dyn = 15;

%% Computing the desired quaternions:

q_d = zeros(4, length(t_list));

for i=1:1:length(t_list)
    q_d(:, i) = eul2quat([yaw_d(i), pitch_d(i), roll_d(i)])';
    q_d(:, i) = [q_d(2, i); q_d(3, i); q_d(4, i); q_d(1, i)];
    %q_d(:, i) = q_d(:, i)/norm(q_d(:, i));

    % if sum(q_d(:, i).^2)~=1.0
    %     sum(q_d(:, i).^2)
    % end
end
q_bar = q_d + (q0 - q_d).*exp(-(t_list - t_list(1))/tau);

for i=1:1:length(t_list)
    normbar = norm(q_bar(:, i));
    q_bar(:, i) = q_bar(:, i)/normbar;
end

%% Computing q_d_dot:

dq_bar = zeros(size(q_bar));
dq_bar(:, 1) = (q_bar(:, 2) - q_bar(:, 1))/dt;

for i=2:1:length(t_list)-1
    dq_bar(:, i) = 0.5*(q_bar(:, i+1) - q_bar(:, i-1))/dt;
end

dq_bar(:, end) = (q_bar(:, end) - q_bar(:, end - 1))/dt;

%% Step 1: Attitude Kinematics:

%Dynamic Inversion:
w_bar = zeros(3, length(t_list));

for i=1:1:length(t_list)
    Gq = G_of_q(q_bar(:, i));
    w_bar(:, i) = 2*Gq*dq_bar(:, i);
end

y0 = [zeros(4, 1); q0-q_d(:,1)];
%Solving for q:
[t_ode, y] = ode45(@(t, y)q_solver(t, y, w_kin), t_list, y0);

%disp(size(y)) = 36001 8
q_err = y(:, 5:8)';
q_err_int = y(:, 1:4)';

w_err = zeros(3, length(t_ode));

for i=1:1:length(t_ode)

    q_bar_t = q_bar(:, i);
    dq_bar_t = dq_bar(:, i);
    
    Gqbar = G_of_q(q_bar_t);
    w_bar_t = 2*Gqbar*dq_bar_t;

    Qwbar = Q_of_w(w_bar_t);
    q_new = q_err(:, i) + q_bar_t;
    normq = norm(q_new);
    q_new = q_new/normq;

    Gq = G_of_q(q_new);

    KI = -2*w_kin*w_kin*Gq;
    KP = -Gq*(4*w_kin*eye(4) + Qwbar);

    Kt = [KI, KP];

    w_err(:, i) = Kt*(y(i, :)');

end

w_calc = w_err + w_bar;
q_calc = q_bar + q_err;

for i=1:1:length(t_list)
    normcalc = norm(q_calc(:, i));
    q_calc(:, i) = q_calc(:, i)/normcalc;
end

%% Step 2: Attitude Dynamics:

w_bar_star = w_calc;
ddt_w_bar_star = zeros(size(w_bar_star));

%Differentiator:
ddt_w_bar_star(:, 1) = (w_bar_star(:,2) - w_bar_star(:,1))/dt;

for i=2:1:length(t_list)-1

    ddt_w_bar_star(:, i) = (w_bar_star(:,i+1) - w_bar_star(:,i-1))/(2*dt);

end

ddt_w_bar_star(:, end) = (w_bar_star(:,end) - w_bar_star(:,end-1))/dt;

%Dynamic Inversion:
Tbar = zeros(size(w_bar_star));

for i=1:1:length(t_list)
    
    Sw = S_of_w(w_bar_star(:, i));
    Tbar(:, i) = Js*ddt_w_bar_star(:, i) + Sw*Js*w_bar_star(:, i);% + Td(:, i);

end

M = M_of_J(Js);

w_guess = [0; 0; 0; w0(1)-w_bar_star(1, 1); w0(2)-w_bar_star(2, 1); w0(3)-w_bar_star(3, 1)];
%w values:
[t_w, y_w] = ode45(@(t, w)w_solver(t, w, w_dyn, M), t_list, w_guess);

%y_w - 181 x 6
w_tilda_star = y_w(:, 4:6)';
w_final = w_bar_star + w_tilda_star;
KI_prime = -w_dyn*w_dyn*Js*M*M;
B_prime = inv(Js); 

T_tilda = zeros(size(w_final));

for i=1:1:length(t_w)

    w_bar_star_t = w_bar_star(:, i);
    w_agm = y_w(i, :);

    w_new = w_final(:, i)+ w_bar_star_t;

    A_prime = -B_prime*(S_of_w(0.5*w_new)*Js + S_of_w(0.5*Js*w_new));
    
    KP_prime = -Js*A_prime - 2*w_dyn*Js*M;

    T_tilda(:, i) = [KI_prime, KP_prime]*w_agm';

end

T_final = T_tilda+Tbar;

%% Step 3: Final value calculations:

q_final = zeros(size(q_bar));
q_final(:, 1) = q0;

for i=2:1:length(t_list)
    
    qi_1 = q_final(:, i-1);

    epsx = W_cross(qi_1(1:3));

    deps = 0.5*(qi_1(4)*eye(3) + epsx)*w_final(:, i-1);
    deta = -0.5*qi_1(1:3)'*w_final(:, i-1);
    
    eps = qi_1(1:3) + dt*deps;
    eta = qi_1(4) + dt*deta;

    qi = [eps; eta];
    normq = norm(qi);
    q_final(:, i) = qi/normq;
end

% for i=2:1:length(t_list)
% 
%     w_t = w_final(:, i-1);
%     Qw = Q_of_w(w_t);
% 
%     dqdt = 0.5*Qw*q_final(:, i-1);
%     q_final(:, i) = q_final(:, i-1) + dqdt*dt;
%     q_final(:, i) = q_final(:, i)/norm(q_final(:, i));
% end
% 
% q_final(:, end) = q_final(:, end-1) + dqdt*dt;
% q_final(:, end) = q_final(:, end)/norm(q_final(:, end));

% [t_final, q_finalT] = ode45(@(t, q)q_final_solver(t, q, t_list, w_final), t_list, q0);
% 
% for i=2:1:length(t_list)
% 
%     w_t = w_final(:, i-1);
%     Qw = Q_of_w(w_t);
% 
%     dq = 0.5*dt*Qw*q_final(:, i-1);
%     dq = [w_t*dt;0];
%     q1 = q_final(:, i-1);
% 
%     % q_fin = quatmultiply([dq(4); dq(1:3)]', [q1(4); q1(1:3)]');
%     % q_final(:, i) = [q_fin(2:4), q_fin(1)]';
%     q_final(:, i) = expm(dt*0.5*Qw)*q_final(:, i-1);
%     q_final(:, i) = q_final(:, i)/norm(q_final(:, i));
% 
% end

% for i=2:1:length(t_list)
% 
%     w_t = w_final(:, i-1);
%     Qw = Q_of_w(w_t);
% 
%     Qw3 = Qw(1:3, :);
% 
%     dq3 = 0.5*Qw3*q_final(:, i-1);
% 
%     qfin = q_final(1:3, i-1) + dq3*dt;
%     qeta = sqrt(1 - qfin'*qfin);
% 
%     q_final(:, i) = [qfin; qeta];
% 
% end

 %% Validation:
eul_found1 = zeros(length(t_list), 3);%quat2eul([0, 1, 2, 3]);
eul_found2 = zeros(length(t_list), 3);

for i=1:1:length(t_list)
    q_temp1 = q_final(:, i);
    q_temp2 = q_calc(:, i);
    q_temp1 = [q_temp1(4), q_temp1(1), q_temp1(2), q_temp1(3)];
    eul_found1(i, :) = quat2eul(q_temp1);
    q_temp2 = [q_temp2(4), q_temp2(1), q_temp2(2), q_temp2(3)];
    eul_found2(i, :) = quat2eul(q_temp2);
end

roll_found1 = eul_found1(:, 3);
pitch_found1 = eul_found1(:, 2); 
yaw_found1 = eul_found1(:, 1);

roll_found2 = eul_found2(:, 3);
pitch_found2 = eul_found2(:, 2); 
yaw_found2 = eul_found2(:, 1);

end

%%

function [T_final, yaw_d, yaw_found1, pitch_d, pitch_found1, roll_d, roll_found1] = standard_TLC(t_ode, tau, c)
%%% Spacecraft model parameters (Wu et al., 2018):

t_list = t_ode;
dt = t_list(2)-t_list(1);

Js = (1+c)*[147, 6.5, 6; 6.5, 158, 5.5; 6, 5.5, 137];
Td = [1+2*sin(0.005*t_list); -1-5*cos(0.005*t_list); 2-4*cos(0.005*t_list)];

w0 = [1; -1; 1]*pi/180;
T0 = [0; 0; 0];

roll0 = 42*pi/180;
pitch0 = -19*pi/180;
yaw0 = -7*pi/180;

q0 = eul2quat([yaw0, pitch0, roll0])';
q0 = [q0(2); q0(3); q0(4); q0(1)];

Tmax = 0.2;

%% Desired trajectory:

if(t_list(end)>1100)

    t1100_ind = find(t_list==1100);
    
    t1 = t_list(1:t1100_ind-1);
    t2 = t_list(t1100_ind:end);
    
    roll_d = [-90*cos(0.025*t1), 0*t2]*pi/180;
    pitch_d = [-30*sin(0.005*t1), 0*t2]*pi/180;
    yaw_d = [75*ones(1, length(t1)), 0*t2]*pi/180;

else

    roll_d = -90*cos(0.025*t_list)*pi/180;
    pitch_d = -30*sin(0.005*t_list)*pi/180;
    yaw_d = 75*ones(1, length(t_list))*pi/180;

end

w_kin = 0.01;
w_dyn = 15;

%% Computing the desired quaternions:

q_d = zeros(4, length(t_list));

for i=1:1:length(t_list)
    q_d(:, i) = eul2quat([yaw_d(i), pitch_d(i), roll_d(i)])';
    q_d(:, i) = [q_d(2, i); q_d(3, i); q_d(4, i); q_d(1, i)];
    %q_d(:, i) = q_d(:, i)/norm(q_d(:, i));

    % if sum(q_d(:, i).^2)~=1.0
    %     sum(q_d(:, i).^2)
    % end
end
q_bar = q_d + (q0 - q_d).*exp(-(t_list - t_list(1))/tau);

% for i=1:1:length(t_list)
%     normbar = norm(q_bar(:, i));
%     q_bar(:, i) = q_bar(:, i)/normbar;
% end

%% Computing q_d_dot:

dq_bar = zeros(size(q_bar));
dq_bar(:, 1) = (q_bar(:, 2) - q_bar(:, 1))/dt;

for i=2:1:length(t_list)-1
    dq_bar(:, i) = 0.5*(q_bar(:, i+1) - q_bar(:, i-1))/dt;
end

dq_bar(:, end) = (q_bar(:, end) - q_bar(:, end - 1))/dt;

%% Step 1: Attitude Kinematics:

%Dynamic Inversion:
w_bar = zeros(3, length(t_list));

for i=1:1:length(t_list)
    Gq = G_of_q(q_bar(:, i));
    w_bar(:, i) = 2*Gq*dq_bar(:, i);
end

y0 = [zeros(4, 1); q0-q_d(:,1)];
%Solving for q:
[t_ode, y] = ode45(@(t, y)q_solver(t, y, w_kin), t_list, y0);

%disp(size(y)) = 36001 8
q_err = y(:, 5:8)';
q_err_int = y(:, 1:4)';

w_err = zeros(3, length(t_ode));

for i=1:1:length(t_ode)

    q_bar_t = q_bar(:, i);
    dq_bar_t = dq_bar(:, i);
    
    Gqbar = G_of_q(q_bar_t);
    w_bar_t = 2*Gqbar*dq_bar_t;

    Qwbar = Q_of_w(w_bar_t);
    % q_new = q_err(:, i) + q_bar_t;
    % normq = norm(q_new);
    % q_new = q_new/normq;

    % Gq = G_of_q(q_new);

    KI = -2*w_kin*w_kin*Gqbar;
    KP = -Gqbar*(4*w_kin*eye(4) + Qwbar);

    Kt = [KI, KP];

    w_err(:, i) = Kt*(y(i, :)');

end

w_calc = w_err + w_bar;
q_calc = q_bar + q_err;

for i=1:1:length(t_list)
    normcalc = norm(q_calc(:, i));
    q_calc(:, i) = q_calc(:, i)/normcalc;
end

%% Step 2: Attitude Dynamics:

w_bar_star = w_calc;
ddt_w_bar_star = zeros(size(w_bar_star));

%Differentiator:
ddt_w_bar_star(:, 1) = (w_bar_star(:,2) - w_bar_star(:,1))/dt;

for i=2:1:length(t_list)-1

    ddt_w_bar_star(:, i) = (w_bar_star(:,i+1) - w_bar_star(:,i-1))/(2*dt);

end

ddt_w_bar_star(:, end) = (w_bar_star(:,end) - w_bar_star(:,end-1))/dt;

%Dynamic Inversion:
Tbar = zeros(size(w_bar_star));

for i=1:1:length(t_list)
    
    Sw = S_of_w(w_bar_star(:, i));
    Tbar(:, i) = Js*ddt_w_bar_star(:, i) + Sw*Js*w_bar_star(:, i);% + Td(:, i);

end

M = M_of_J(Js);

w_guess = [0; 0; 0; w0(1)-w_bar_star(1, 1); w0(2)-w_bar_star(2, 1); w0(3)-w_bar_star(3, 1)];
%w values:
[t_w, y_w] = ode45(@(t, w)w_solver(t, w, w_dyn, M), t_list, w_guess);

%y_w - 181 x 6
w_tilda_star = y_w(:, 4:6)';
w_final = w_bar_star + w_tilda_star;
KI_prime = -w_dyn*w_dyn*Js*M*M;
B_prime = inv(Js); 

T_tilda = zeros(size(w_final));

for i=1:1:length(t_w)

    w_bar_star_t = w_bar_star(:, i);
    w_agm = y_w(i, :);

    %w_new = w_final(:, i)+ w_bar_star_t;

    A_prime = -B_prime*(S_of_w(w_bar_star_t)*Js + S_of_w(Js*w_bar_star_t));
    
    KP_prime = -Js*A_prime - 2*w_dyn*Js*M;

    T_tilda(:, i) = [KI_prime, KP_prime]*w_agm';

end

T_final = T_tilda+Tbar;

%% Step 3: Final value calculations:

q_final = zeros(size(q_bar));
q_final(:, 1) = q0;

for i=2:1:length(t_list)
    
    qi_1 = q_final(:, i-1);

    epsx = W_cross(qi_1(1:3));

    deps = 0.5*(qi_1(4)*eye(3) + epsx)*w_final(:, i-1);
    deta = -0.5*qi_1(1:3)'*w_final(:, i-1);
    
    eps = qi_1(1:3) + dt*deps;
    eta = qi_1(4) + dt*deta;

    qi = [eps; eta];
    normq = norm(qi);
    q_final(:, i) = qi/normq;
end

% for i=2:1:length(t_list)
% 
%     w_t = w_final(:, i-1);
%     Qw = Q_of_w(w_t);
% 
%     dqdt = 0.5*Qw*q_final(:, i-1);
%     q_final(:, i) = q_final(:, i-1) + dqdt*dt;
%     q_final(:, i) = q_final(:, i)/norm(q_final(:, i));
% end
% 
% q_final(:, end) = q_final(:, end-1) + dqdt*dt;
% q_final(:, end) = q_final(:, end)/norm(q_final(:, end));

% [t_final, q_finalT] = ode45(@(t, q)q_final_solver(t, q, t_list, w_final), t_list, q0);
% 
% for i=2:1:length(t_list)
% 
%     w_t = w_final(:, i-1);
%     Qw = Q_of_w(w_t);
% 
%     dq = 0.5*dt*Qw*q_final(:, i-1);
%     dq = [w_t*dt;0];
%     q1 = q_final(:, i-1);
% 
%     % q_fin = quatmultiply([dq(4); dq(1:3)]', [q1(4); q1(1:3)]');
%     % q_final(:, i) = [q_fin(2:4), q_fin(1)]';
%     q_final(:, i) = expm(dt*0.5*Qw)*q_final(:, i-1);
%     q_final(:, i) = q_final(:, i)/norm(q_final(:, i));
% 
% end

% for i=2:1:length(t_list)
% 
%     w_t = w_final(:, i-1);
%     Qw = Q_of_w(w_t);
% 
%     Qw3 = Qw(1:3, :);
% 
%     dq3 = 0.5*Qw3*q_final(:, i-1);
% 
%     qfin = q_final(1:3, i-1) + dq3*dt;
%     qeta = sqrt(1 - qfin'*qfin);
% 
%     q_final(:, i) = [qfin; qeta];
% 
% end

 %% Validation:
eul_found1 = zeros(length(t_list), 3);%quat2eul([0, 1, 2, 3]);
eul_found2 = zeros(length(t_list), 3);

for i=1:1:length(t_list)
    q_temp1 = q_final(:, i);
    q_temp2 = q_calc(:, i);
    q_temp1 = [q_temp1(4), q_temp1(1), q_temp1(2), q_temp1(3)];
    eul_found1(i, :) = quat2eul(q_temp1);
    q_temp2 = [q_temp2(4), q_temp2(1), q_temp2(2), q_temp2(3)];
    eul_found2(i, :) = quat2eul(q_temp2);
end

roll_found1 = eul_found1(:, 3);
pitch_found1 = eul_found1(:, 2); 
yaw_found1 = eul_found1(:, 1);

roll_found2 = eul_found2(:, 3);
pitch_found2 = eul_found2(:, 2); 
yaw_found2 = eul_found2(:, 1);

end

%% Function wx:

function wx = W_cross(w)
    
    wx = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];

end
%% Function G(q):

function G = G_of_q(q)
    
    G = [q(4), q(3), -q(2), -q(1); -q(3), q(4), q(1), -q(2); q(2), -q(1), q(4), -q(3)];

end

%% Function Q(w):

function Q = Q_of_w(w)
    
    Q = [0, w(3), -w(2), w(1); -w(3), 0, w(1), w(2); w(2), -w(1), 0, w(3); -w(1), -w(2), -w(3), 0];

end

%% Function S(w):

function S = S_of_w(w)
    
    S = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];

end

%% Function M:

function M = M_of_J(J)
    
    MatJ = diag([sqrt(J(1,1)), sqrt(J(2,2)), sqrt(J(3,3))]);
    M = inv(MatJ); %diag([1/sqrt(J(1,1)), 1/sqrt(J(2,2)), 1/sqrt(J(3,3))]);

end

%% Function q error:

function dqdt = q_solver(t, q, w_kin)
    
  dqdt = [zeros(4), eye(4); -w_kin*w_kin*eye(4), -2*w_kin*eye(4)]*q;

end

%% Function w error:

function dwdt = w_solver(t, w, w_dyn, M)

    dwdt = [zeros(3), eye(3); -w_dyn*w_dyn*M*M, -2*w_dyn*M]*w;

end

%% Function q final:

function dqdt = q_final_solver(t, q, t_list, w_list)

    w_t1 = interp1(t_list, w_list(1, :), t);
    w_t2 = interp1(t_list, w_list(2, :), t);
    w_t3 = interp1(t_list, w_list(3, :), t);

    w_t = [w_t1; w_t2; w_t3];

    Qw = Q_of_w(w_t);

    dqdt = 0.5*Qw*q;

end
