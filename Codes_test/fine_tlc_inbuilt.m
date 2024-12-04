%% Fine TLC (Inbuilt)

clc; clear; close all;

%% Spacecraft model parameters (Wu et al., 2018):

t_list = 0:0.5:1100;
dt = t_list(end) - t_list(end-1); 

Js = [147, 6.5, 6; 6.5, 158, 5.5; 6, 5.5, 137];
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
tao = 1;

%% Computing the desired quaternions:

q_d = zeros(4, length(t_list));

for i=1:1:length(t_list)
    q_d(:, i) = eul2quat([yaw_d(i), pitch_d(i), roll_d(i)])';
    q_d(:, i) = [q_d(2, i); q_d(3, i); q_d(4, i); q_d(1, i)];
end
q_bar = q_d + (q0 - q_d).*exp(-(t_list - t_list(1))/tao);

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

y0 = [q0, zeros(4, 1)];
%Solving for q:
[t_ode, y] = ode45(@(t, y)q_solver(t, y, w_kin, q0, t_list(1), tao), t_list, y0);

%disp(size(y)) = 36001 8
q_err = y(:, 1:4)';
q_err_int = y(:, 5:8)';

w_err = zeros(3, length(t_ode));

for i=1:1:length(t_ode)

    q_bar_t = q_bar(:, i);
    dq_bar_t = dq_bar(:, i);
    
    Gqbar = G_of_q(q_bar_t);
    w_bar_t = 2*Gqbar*dq_bar_t;

    Qwbar = Q_of_w(w_bar_t);
    Gq = G_of_q(q_err(:, i) + q_bar_t);

    KI = -2*w_kin*w_kin*Gq;
    KP = -Gq*(4*w_kin*eye(4) + Qwbar);

    w_err(:, i) = KI*q_err_int(:, i) + KP*q_err(:, i);

end

w_calc = w_err + w_bar;
q_calc = q_bar + q_err;

%% Step 2: Attitude Dynamics:

w_bar_star = w_calc;
ddt_w_bar_star = zeros(size(w_bar_star));

%Differentiator:
syms t s

wd = 0.25;
xi = 0.707;

D2 = wd*wd*s/(s^2 + 2*wd*xi*s + wd*wd);

n = 15;%(floor((1/(dt))));

[Snx1, ak_list1, bk_list1] = fourier(w_bar_star(1, :), t_list, n);
[Snx2, ak_list2, bk_list2] = fourier(w_bar_star(2, :), t_list, n);
[Snx3, ak_list3, bk_list3] = fourier(w_bar_star(3, :), t_list, n);

figure(100)
plot(t_list, w_bar_star(1, :))
hold on
plot(t_list, Snx1)
plot(t_list, w_bar_star(2, :))
plot(t_list, Snx2)
plot(t_list, w_bar_star(3, :))
plot(t_list, Snx3)

L = t_list(end) - t_list(1);
ft1 = fourier_series(t, ak_list1, bk_list1, L);
ft2 = fourier_series(t, ak_list2, bk_list2, L);
ft3 = fourier_series(t, ak_list3, bk_list3, L);
w_bar_star_timefn = [ft1; ft2; ft3];

lft1 = laplace(ft1);
lft2 = laplace(ft2);
lft3 = laplace(ft3);

laplace_w = [lft1; lft2; lft3];
ddt_w_bar_star_freqdomain = D2*laplace_w;
ddt_w_bar_star_timedomain = ilaplace(ddt_w_bar_star_freqdomain);

for i=1:1:length(t_list)

    ddt_w_bar_star(:, i) = double(subs(ddt_w_bar_star_timedomain, t, t_list(i)));

end

%Dynamic Inversion:
Tbar = zeros(size(w_bar_star));

for i=1:1:length(t_list)
    
    Sw = S_of_w(w_bar_star(:, i));
    Tbar(:, i) = Js*ddt_w_bar_star(:, i) + Sw*Js*w_bar_star(:, i);

end

w_guess = [0; 0; 0; w0(1)-w_bar_star(1, 1); w0(2)-w_bar_star(2, 1); w0(3)-w_bar_star(3, 1)];
%w values:
[t_w, y_w] = ode45(@(time, w)w_solver(time, w, w_dyn, w_bar_star_timefn, Js, t), t_list, w_guess);

%y_w - 181 x 6
w_tilda_star = y_w(:, 4:6);
w_final = w_bar_star + w_tilda_star';
KI_prime = -w_dyn*w_dyn*eye(3);
B_prime = diag([1/Js(1,1), 1/Js(2,2), 1/Js(3,3)]); 

T_tilda = zeros(size(w_final));

for i=1:1:length(t_w)

    w_bar_star_t = w_bar_star(:, i);
    w_agm = y_w(i, :);

    w_new = w_final(:, i)+ w_bar_star_t;

    A_prime = -B_prime*(S_of_w(0.5*w_new)*Js + S_of_w(0.5*Js*w_new));
    KP_prime = -Js*A_prime - 2*w_dyn*sqrt(Js);

    T_tilda(:, i) = [KI_prime, KP_prime]*w_agm';

end

T_final = T_tilda+Tbar;

%% Step 3: Final value calculations:

%%%
[Sf1, ak_listf1, bk_listf1] = fourier(w_final(1, :), t_list, n);
[Sf2, ak_listf2, bk_listf2] = fourier(w_final(2, :), t_list, n);
[Sf3, ak_listf3, bk_listf3] = fourier(w_final(3, :), t_list, n);

figure(101)
plot(t_list, w_final(1, :))
hold on
plot(t_list, Sf1)
plot(t_list, w_final(2, :))
plot(t_list, Sf2)
plot(t_list, w_final(3, :))
plot(t_list, Sf3)

L = t_list(end) - t_list(1);
wft1 = fourier_series(t, ak_listf1, bk_listf1, L);
wft2 = fourier_series(t, ak_listf2, bk_listf2, L);
wft3 = fourier_series(t, ak_listf3, bk_listf3, L);
w_final_time = [wft1; wft2; wft3];
%%%
 
[t_final, q_final] = ode45(@(time, q)q_final_solver(time, q, w_final_time, t), t_list, q0);


%% Validation:
eul_found = zeros(length(t_list), 3);%quat2eul([0, 1, 2, 3]);

for i=1:1:length(t_list)
    q_temp = q_final(i, :);
    q_temp = [q_temp(4), q_temp(1), q_temp(2), q_temp(3)];
    eul_found(i, :) = quat2eul(q_temp);
end

roll_found = eul_found(:, 3);
pitch_found = eul_found(:, 2); 
yaw_found = eul_found(:, 1);


%% Plotting Block:

%Roll:
figure(1)
plot(t_ode, roll_d*180/pi, 'LineWidth', 2, 'LineStyle','--', 'Color', 'g', 'DisplayName','Desired')
hold on
plot(t_ode, roll_found*180/pi, 'LineWidth', 2, 'LineStyle','-', 'Color', 'b', 'DisplayName','Fine TLC')
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('Roll angle')
xlabel("Time (s) -->")
ylabel("Roll (deg) -->")
%axis padded

%Pitch:
figure(2)
plot(t_ode, pitch_d*180/pi, 'LineWidth', 2, 'LineStyle','--', 'Color', 'g', 'DisplayName','Desired')
hold on
plot(t_ode, pitch_found*180/pi, 'LineWidth', 2, 'LineStyle','-', 'Color', 'b', 'DisplayName','Fine TLC')
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('Pitch angle')
xlabel("Time (s) -->")
ylabel("Pitch (deg) -->")

%Yaw:
figure(3)
plot(t_ode, yaw_d*180/pi, 'LineWidth', 2, 'LineStyle','--', 'Color', 'g', 'DisplayName','Desired')
hold on
plot(t_ode, yaw_found*180/pi, 'LineWidth', 2, 'LineStyle','-', 'Color', 'b', 'DisplayName','Fine TLC')
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('Yaw angle')
xlabel("Time (s) -->")
ylabel("Yaw (deg) -->")

%Omega:
figure(4)
plot(t_ode, abs(w_calc(1, :)), 'LineWidth', 2, 'LineStyle','--', 'Color', 'g', 'DisplayName','Desired')
hold on
%plot(t, yaw_found*180/pi, 'LineWidth', 2, 'LineStyle','-', 'Color', 'b', 'DisplayName','Fine TLC')
legend("show")
ax = gca;
ax.FontSize = 16;
grid on
title('Yaw angle')
xlabel("Time (s) -->")
ylabel("Yaw (deg) -->")



%% Function G(q):

function G = G_of_q(q)
    
    G = [q(4), q(3), -q(2), -q(1); -q(3), q(4), q(1), -q(2); q(2), -q(1), q(4), -q(3)];
    %G = [-q(2), q(1), q(4), -q(3); -q(3), -q(4), q(1), q(2); -q(4), q(3), -q(2), q(1)];

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
    
    M = diag([1/sqrt(J(1,1)), 1/sqrt(J(2,2)), 1/sqrt(J(3,3))]);

end

%% Function q error:

function dqdt = q_solver(t, q, w_kin, q0, t0, tao)

    dqdt = zeros(size(q));
    dqdt(5:8) = q(1:4);
    dt = 0.001;
 
    if(t<1100)

        roll_d = -90*cos(0.025*t)*pi/180;
        pitch_d = -30*sin(0.005*t)*pi/180;
        yaw_d = 75*pi/180;

    else

        roll_d = 0;
        pitch_d = 0;
        yaw_d = 0;

    end

    if(t<1100)

        roll_ddt = -90*cos(0.025*(t+dt))*pi/180;
        pitch_ddt = -30*sin(0.005*(t+dt))*pi/180;
        yaw_ddt = 75*pi/180;

    else

        roll_ddt = 0;
        pitch_ddt = 0;
        yaw_ddt = 0;

    end

    q_bar_t = eul2quat([yaw_d, pitch_d, roll_d])';
    q_bar_t = [q_bar_t(2); q_bar_t(3); q_bar_t(4); q_bar_t(1)];
    q_bar_t = q_bar_t + (q0 - q_bar_t).*exp(-(t - t0)/tao);

    q_bar_tdt = eul2quat([yaw_ddt, pitch_ddt, roll_ddt])';
    q_bar_tdt = [q_bar_tdt(2); q_bar_tdt(3); q_bar_tdt(4); q_bar_tdt(1)];
    q_bar_tdt = q_bar_tdt + (q0 - q_bar_tdt).*exp(-(t+dt - t0)/tao);

    dq_bar_t = (q_bar_tdt - q_bar_t)/dt;

    Gqbar = G_of_q(q_bar_t);
    w_bar_t = 2*Gqbar*dq_bar_t;

    Qwbar = Q_of_w(w_bar_t);
    Gq = G_of_q(q(1:4) + q_bar_t);

    KI = -2*w_kin*w_kin*Gq;
    KP = -Gq*(4*w_kin*eye(4) + Qwbar);

    K = [KI, KP];

    dqdt(1:4) = 0.5*Qwbar*(q(1:4)) + 0.5*Gq'*K*[q(5:8); q(1:4)];

end

%% Function w error:

function dwdt = w_solver(time, w, w_dyn, w_bar_star, Js, t)

    %dwdt(1:3) = q(4:6);
    
    w_bar_star_t = double(subs(w_bar_star, t, time));

    w_new = w(4:6)+2*w_bar_star_t;

    B_prime = inv(Js); 
    A_prime = -B_prime*(S_of_w(0.5*w_new)*Js + S_of_w(0.5*Js*w_new));

    KI_prime = -w_dyn*w_dyn*eye(3);
    KP_prime = -Js*A_prime - 2*w_dyn*sqrt(Js);

    Mat1 = [zeros(3), eye(3); zeros(3), A_prime];
    Mat2 = [zeros(3); B_prime]*[KI_prime, KP_prime];

    dwdt = (Mat1 + Mat2)*w;

end

%% Function w error:

function dqdt = q_final_solver(time, q, w_final_time, t)

    w_t = double(subs(w_final_time, t, time));

    Qw = Q_of_w(w_t);

    dqdt = 0.5*Qw*q;

end

%% Fourier Series:

function [Snx, ak_list, bk_list] = fourier(fx, x, n)

    L = x(end) - x(1);
    a0 = trapz(x, fx);
    a0 = 2*a0/L;

    ak_list = zeros(1, n+1);
    ak_list(1) = a0;
    bk_list = zeros(1, n);

    for k=1:1:n
        
        ak = trapz(x, fx.*cos(2*pi*k*x/L));
        ak = 2*ak/L;
        ak_list(k) = ak;

        bk = trapz(x, fx.*sin(2*pi*k*x/L));
        bk = 2*bk/L;
        bk_list(k) = bk;

    end

    Snx = a0/2;

    for k=1:1:n

        Snx = Snx + ak_list(k+1)*cos(2*pi*k*x/L) + bk_list(k)*sin(2*pi*k*x/L);

    end

end

%% Fourier function:

function ft = fourier_series(t, ak_list, bk_list, L)

    ft = ak_list(1)/2;

    for i=1:1:length(bk_list)

        ft = ft + ak_list(i)*cos(2*pi*i*t/L) + bk_list(i)*sin(2*pi*i*t/L);

    end

end