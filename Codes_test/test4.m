%% Fine TLC (Inbuilt)

clc; clear; close all;

%% Spacecraft model parameters (Wu et al., 2018):

n=100;
l = 2*n + 2;

t_list = linspace(0,2*pi, l);
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
t_mod = t_list - t_list(1);

q_bar = q_d + (q0 - q_d).*exp(-(t_mod)/tao);

figure(1)
plot(t_list, q_d + (q0 - q_d).*exp(-t_mod))%exp(-(t_mod)/tao))
hold on
% plot(t_list, q_bar(2, :), 'g')
% plot(t_list, q_bar(3, :))
% plot(t_list, q_bar(4, :))

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

Px1 = Lagrange(50, t_list, 1, l, 0, 2*pi, w_bar(1, :)');

[Snx1, ak_list1, bk_list1] = fourier(w_bar(1, :), t_list, n);
[Snx2, ak_list2, bk_list2] = fourier(w_bar(2, :), t_list, n);
[Snx3, ak_list3, bk_list3] = fourier(w_bar(3, :), t_list, n);

figure(100)
plot(t_list, w_bar(1, :))
hold on
plot(t_list, Snx1, 'g')
% plot(t_list, w_bar(2, :))
% plot(t_list, Snx2)
% plot(t_list, w_bar(3,:))
% plot(t_list, Snx3)

%% Function G(q):

function G = G_of_q(q)
    
    G = [q(4), q(3), -q(2), -q(1); -q(3), q(4), q(1), -q(2); q(2), -q(1), q(4), -q(3)];
    %G = [-q(2), q(1), q(4), -q(3); -q(3), -q(4), q(1), q(2); -q(4), q(3), -q(2), q(1)];

end

%% Fourier Series:

function [Snx, ak_list, bk_list] = fourier(fx, x, n)

    L = 2*pi;%x(end) - x(1);
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

%% Returns the error array:

function Px = Lagrange(n, x_vals, nodes_fn, N, lb, ub, fx)
    
    if nodes_fn == 1
        x_nodes = linspace(lb, ub, n+1);
    end

    if nodes_fn == 2
        x_nodes = chebyshev(n, lb, ub);
    end
    
    Lk = zeros(n+1, N);
    
    for k=1:1:n+1
        den = prod(x_nodes(k) - [x_nodes(1:k-1), x_nodes(k+1:end)]);
    
        for i = 1:1:N
            num = (prod(x_vals(i) - [x_nodes(1:k-1) x_nodes(k+1:n+1)]));
            Lk(k, i) = num/den;
        end
        
    end
    
    Px = zeros(1, N);
    
    for i=1:1:N
       % x = x_vals(i);
        Px(1, i) = sum(Lk(:, i).*fx);
    end

end 
