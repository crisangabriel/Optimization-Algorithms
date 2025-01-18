
Q = [1, -2; -3, 7]; % Corrected symmetric matrix
Q = 0.5 * (Q + Q')
c = [2; 4];
A = [1, -0.28; 1, 0.33; -1, 0.6; -1, -0.4];
b = [0.57; 3.66; 3.8; 0.8];
%x0 = [2; 5]; % Initial feasible point
%W0 = [2;5]; % Initial working set
max_iter = 1000;
toll = 1e-6;

[x, info] = ASQP_withfindingconstraints(Q, c, A, b, max_iter, toll);

