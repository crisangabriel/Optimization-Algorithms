
% Example usage:
Q = [8, 5.5; 5.5, 13];
c = [-7; 2];
A = [1, 10; 1, -1; -1, 3; -1, -3];
b = [25; 3; 4; -6];
x0 = [25/13; 29/13];
W0 = [1; 3];
max_iter = 1000;
toll = 1e-6;

[x, info] = ASQPv2(Q, c, A, b, x0, W0, max_iter, toll);