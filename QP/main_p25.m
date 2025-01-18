clc;
clear all;
%% Example usage:
Q = [1, -2; -3, 7]; % Provided matrix
Q = 0.5 * (Q + Q'); % Make it symmetric
c = [2; 4];
A = [1, -0.28; 1, 0.33; -1, 0.6; -1, -0.4];
b = [0.57; 3.66; 3.8; 0.8];
x0 = [-2; 3]; % Provided initial feasible point
max_iter = 1000;
toll = 1e-6;

[x, info] = ASQP_prob25(Q, c, A, b, x0, max_iter, toll);