function [x_opt, f_opt] = rosenbrock_method(f, x0, d, s, alpha, beta, epsilon, max_iterations)
    % Initialize variables
    n = length(x0);
    k = 1;
    x = x0;
    stop_criterion = false;
    
    % Initialize plot
    figure;
    hold on;
    plot(x(1), x(2), 'ro', 'MarkerSize', 8, 'DisplayName', 'Start Point');
    xlabel('x1');
    ylabel('x2');
    title('Progress of Rosenbrock Method');
  
    
    % Store points for plotting
    x_points = x';
    
    % Main loop
    while ~stop_criterion && k <= max_iterations
        % Initialize vectors
        c = zeros(n, 1);
        success = zeros(n, 1);
        fail = zeros(n, 1);
        oscillation = false;

        while ~oscillation
            for i = 1:n
                fprintf('Iteration %d, Direction %d\n', k, i);
                fprintf('Current point: [%s]\n', num2str(x'));
                fprintf('Step size: %f\n', s(i));

                x_new = x + s(i) * d(:,i);
                f_x_new = f(x_new);
                f_x = f(x);

                if f_x_new <= f_x
                    x = x_new;
                    k = k + 1;
                    success(i) = 1;
                    c(i) = c(i) + s(i);
                    s(i) = s(i) * alpha;
                    fprintf('Success in direction %d\n', i);
                else
                    fail(i) = 1;
                    s(i) = s(i) * beta;
                    fprintf('Failure in direction %d\n', i);
                    fprintf('Condition: f(x + s * d) > f(x)\n');
                    fprintf('f(x + s * d) = %f, f(x) = %f\n', f_x_new, f_x);
                end

                fprintf('New point: [%s]\n', num2str(x'));
                fprintf('Objective function value: %f\n', f(x));

                % Store the new point and plot the line
                x_points = [x_points; x'];
                plot(x_points(:,1), x_points(:,2), 'b-');
                plot(x(1), x(2), 'b.', 'MarkerSize', 8);
                drawnow;
            end

            if all(success == 1) && all(fail == 1)
                oscillation = true;
                fprintf('Oscillation detected. Stopping.\n');
                fprintf('Values of s after oscillation: [%s]\n', num2str(s'));
                fprintf('Values of c after oscillation: [%s]\n', num2str(c'));
                fprintf('Last point calculated: [%s]\n', num2str(x'));
                break;
            end
        end

        if oscillation
            % Compute new directions using the Gram-Schmidt process
            a = d * c; % This results in a vector of length n
            b = zeros(size(d)); % Initialize b with the same size as d

            fprintf('Computing new directions using the Gram-Schmidt procedure:\n');

            for i = 1:n
                b(:,i) = a; 
                for j = 1:i-1
                    % Add check to prevent division by zero
                    if norm(d(:,j)) == 0
                        continue;
                    end
                    b(:,i) = b(:,i) - (a' * d(:,j)) / (d(:,j)' * d(:,j)) * d(:,j);
                end

                % Add check to prevent division by zero
                if norm(b(:,i)) == 0
                    continue;
                end
                d(:,i) = b(:,i) / norm(b(:,i));

                % Print computation details
                fprintf('a_%d = [%s]\n', i, num2str(a'));
                fprintf('b_%d = [%s]\n', i, num2str(b(:,i)'));
                fprintf('New direction d_%d = [%s]\n', i, num2str(d(:,i)'));
            end
        end

        % Check stop criterion
        if norm(s) < epsilon
            stop_criterion = true;
        end
    end

    x_opt = x;
    f_opt = f(x);

    if k > max_iterations
        fprintf('Maximum number of iterations reached: %d\n', max_iterations);
    end

    % Final plot adjustments
    plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 8, 'DisplayName', 'Optimal Point');
    hold off;
end
