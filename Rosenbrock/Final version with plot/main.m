% Objective function
f = @(x) (x(1)^2 + 2*x(1)*x(2)+10*x(2)^2+45*x(2)^6);

% Initial point
x0 = [4, 4]';

% Initial directions
d = eye(length(x0));

% Step lengths
s = [1, 2]';

% Parameters
alpha = 3; 
beta = -0.8;
epsilon = 1e-3; % Tolerance
max_iterations = 20;

[x_opt, f_opt] = rosenbrock_method(f, x0, d, s, alpha, beta, epsilon, max_iterations);

% Display the results
disp('Optimal point:')
disp(x_opt)
disp('Optimal function value:')
disp(f_opt)
