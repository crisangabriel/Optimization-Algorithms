clear all
clc
%%
Q = [2 , 0; 0, 2]; % Corrected symmetric matrix
c = [-4; 6];
A = [1, 1; 1, -1; -1, 1; -1, -1];
b = [2; 2; 2; 2];
x0 = [2; 0]; % Initial feasible point
W0 = [1,2]; % Initial working set
max_iter = 1000;
toll = 1e-6;

[x, info] = ASQP_prob27(Q, c, A, b, x0, W0, max_iter, toll);