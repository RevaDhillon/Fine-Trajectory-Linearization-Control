clc; clear; close all;

% syms s t
% 
% wd = 0.25;
% xi = 0.707;
% 
% D2 = wd*wd*s/(s^2 + 2*wd*xi*s + wd*wd);
% ft = ilaplace(D2);
% 
% double(subs(ft, t, 0))

q = [1/sqrt(2), 0, 0, 1/sqrt(2)];
eul = quat2eul(q)

eul_test = [0, 0, pi/2];
q_test = eul2quat(eul_test)