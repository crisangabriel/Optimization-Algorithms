% Define the parameters for the quadratic programming problem

% Quadratic term matrix Q
Q = [1, 4; 1, 15];

% Linear term vector c
c = [-5; -1];

% Constraint matrix A
A = [1, -0.4;
     1, 0.25;
     -1, 1;
     -1, -0.6];

% Constraint bounds b
b = [1.4; 4; 6; -0.4];

% Initial feasible point x0 (must satisfy Ax <= b)
x0 = [0; 0];

% Initial working set W0 (indices of active constraints, starting with none)
W0 = [];

% Run the active set method
active_set_method(Q, c, A, b, x0, W0);
