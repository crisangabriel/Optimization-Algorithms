function active_set_method(Q, c, A, b, x0, W0)
    % Input parameters
    % Q: Quadratic term matrix
    % c: Linear term vector
    % A: Constraint matrix
    % b: Constraint bounds
    % x0: Initial feasible point
    % W0: Initial working set (indices of active constraints)

    % Initial setup
    xk = x0;
    Wk = W0;
    g0 = Q * x0 + c;
    A_wk = A(Wk, :);

    % Solve the linear system to find initial search direction
    [dk, lambda] = solve_linear_system(Q, A_wk, g0);
    
    % To store all points for plotting
    points = xk';
    
    % Iteration counter
    k = 0;
    
    max_iterations = 1000; % Safety check to avoid infinite loops

    while k < max_iterations
        % Check stopping criteria
        if all(lambda >= 0) || norm(dk) < 1e-6
            fprintf('Optimal solution found at iteration %d\n', k);
            disp('Optimal point:');
            disp(xk);
            break;
        end

        % Determine step length alpha_k
        if any(lambda < 0)
            [~, j] = min(lambda);
            fprintf('Lambda: ');
            disp(lambda);
            fprintf('Removing constraint %d from working set\n', Wk(j));
            Wk(j) = [];
        else
            alpha_k = compute_step_length(A, b, xk, dk, Wk);
            fprintf('Step length alpha_k: %.6f\n', alpha_k);
            xk = xk + alpha_k * dk;
            points = [points; xk'];

            if alpha_k < 1
                ib = find_blocking_constraint(A, b, xk, dk, Wk);
                Wk = [Wk; ib];
                fprintf('Adding constraint %d to working set\n', ib);
            end
        end
        
        % Update iteration counter
        k = k + 1;

        % Compute gradient at new point
        gk = Q * xk + c;

        % Solve the linear system to find new search direction
        A_wk = A(Wk, :);
        [dk, lambda] = solve_linear_system(Q, A_wk, gk);

        % Diagnostic printout for debugging
        fprintf('Iteration %d\n', k);
        disp('Current point:');
        disp(xk);
        disp('Current search direction:');
        disp(dk);
        disp('Current lambda:');
        disp(lambda);
    end

    if k == max_iterations
        fprintf('Reached maximum number of iterations without convergence.\n');
    end

    % Plot the points
    figure;
    plot(points(:, 1), points(:, 2), '-o');
    xlabel('x1');
    ylabel('x2');
    title('Active Set Method Iterations');
    grid on;
end

function [dk, lambda] = solve_linear_system(Q, A_wk, g)
    % Solve the KKT system
    n = size(Q, 1);
    m = size(A_wk, 1);
    KKT = [Q, A_wk'; A_wk, zeros(m, m)];
    rhs = -[g; zeros(m, 1)];
    sol = KKT \ rhs;
    dk = sol(1:n);
    lambda = sol(n+1:end);
end

function alpha_k = compute_step_length(A, b, xk, dk, Wk)
    alpha_k = 1;
    for i = 1:size(A, 1)
        if ~ismember(i, Wk)
            ai = A(i, :);
            if ai * dk > 0
                alpha_i = (b(i) - ai * xk) / (ai * dk);
                alpha_k = min(alpha_k, alpha_i);
            end
        end
    end
end

function ib = find_blocking_constraint(A, b, xk, dk, Wk)
    alpha_k = 1;
    ib = -1;
    for i = 1:size(A, 1)
        if ~ismember(i, Wk)
            ai = A(i, :);
            if ai * dk > 0
                alpha_i = (b(i) - ai * xk) / (ai * dk);
                if alpha_i < alpha_k
                    alpha_k = alpha_i;
                    ib = i;
                end
            end
        end
    end
end