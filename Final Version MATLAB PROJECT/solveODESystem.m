function [t, y] = solveODESystem(f, tspan, y0, N)
    % solveODESystem: Solves a system of ordinary differential equations using RK4 method.
    % Inputs:
    %   - f: Function handle representing the system of ODEs, f(t, y).
    %   - tspan: Time span [t0, tf].
    %   - y0: Initial conditions at t0.
    %   - N: Number of time steps.
    % Outputs:
    %   - t: Array of time points.
    %   - y: Array of solution values at corresponding time points.

    % Initialize parameters
    t0 = tspan(1);
    tf = tspan(2);
    h = (tf - t0) / N;

    % Initialize arrays
    t = linspace(t0, tf, N+1);
    y = zeros(length(y0), N+1);
    y(:, 1) = y0;

    % Perform RK4 method
    for i = 1:N
        k1 = h * f(t(i), y(:, i));
        k2 = h * f(t(i) + 0.5 * h, y(:, i) + 0.5 * k1);
        k3 = h * f(t(i) + 0.5 * h, y(:, i) + 0.5 * k2);
        k4 = h * f(t(i) + h, y(:, i) + k3);

        y(:, i+1) = y(:, i) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end

