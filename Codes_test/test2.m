
clc; close all; clear;

t_list = 0:0.0001:0.02;
w_bar_star = load("w_bar_star.mat");
w_bar_star = w_bar_star.w_bar_star;

[Snx1, ak_list1, bk_list1] = fourier(w_bar_star(1, :), t_list, 1000);
[Snx2, ak_list2, bk_list2] = fourier(w_bar_star(2, :), t_list, 1000);
[Snx3, ak_list3, bk_list3] = fourier(w_bar_star(3, :), t_list, 1000);

figure(100)
plot(t_list, w_bar_star(1, :), 'k')
hold on
plot(t_list, Snx1, 'k')
plot(t_list, w_bar_star(2, :), 'b')
plot(t_list, Snx2, 'b')
plot(t_list, w_bar_star(3, :), 'r')
plot(t_list, Snx3, 'r')

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