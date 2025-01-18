% Objective function
f = @(x) (x(1)^2 - 2*x(1) + 2*x(2) + x(2)^2);

% Initial point
x0 = [0, 1]';

% Initial directions
d = eye(length(x0));

% Step lengths
s = [0.5, -0.5]';

% Parameters
alpha = 2; 
beta = -0.8;
epsilon = 1e-3; % Tolerance
max_iterations = 20;

[x_opt, f_opt] = rosenbrock_method(f, x0, d, s, alpha, beta, epsilon, max_iterations);

% Display the results
disp('Optimal point:')
disp(x_opt)
disp('Optimal function value:')
disp(f_opt)

%% Oscillation detected. Stopping.
Values of s after oscillation: [-1.6         1.6]
Values of c after oscillation: [1.5        -1.5]
Last point calculated: [1.5        -0.5]
Computing new directions using the Gram-Schmidt procedure:
a_1 = [1.5        -1.5]
b_1 = [1.5        -1.5]
New direction d_1 = [0.70711    -0.70711]