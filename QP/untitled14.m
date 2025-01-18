% Define the parameters for the quadratic programming problem

% Quadratic term matrix Q
Q = [8, 5.5; 5.5, 13];

% Linear term vector c
c = [-7; 2];

% Constraint matrix A
A = [1, 10;
     1, -1;
     -1, 3;
     -1, -3];

% Constraint bounds b
b = [25; 3; 4; -6];

% Initial feasible point x0 (intersection of C1 and C3)
x0 = [25/13; 29/13];

% Initial working set W0 (indices of active constraints C1 and C3)
W0 = [1; 3];

% Run the active set method
active_set_method(Q, c, A, b, x0, W0);