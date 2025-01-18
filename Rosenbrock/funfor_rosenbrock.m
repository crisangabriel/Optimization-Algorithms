clear all;
clc;
%%
% Objective function
f = @(x) (x(1)^2 - 2*x(1) + 2*x(2) + x(2)^2);

% Initial point
x0 = [2, 2]';

% Initial directions
d = eye(length(x0));

% Step lengths
s = [0.5, -0.5]';

% Parameters
alpha = 2; 
beta = -0.8;
epsilon = 1e-3; % Tolerance
max_iterations = 20;% Number if iterations

[x_opt, f_opt] = rosenbrock_methodv3osc_good(f, x0, d, s, alpha, beta, epsilon,max_iterations);

% Display the results
disp('Optimal point:')
disp(x_opt)
disp('Optimal function value:')
disp(f_opt)