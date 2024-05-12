function PDE_Solver_GUI()
    % Create a figure
    fig = figure('Position', [200, 200, 600, 400], 'Name', 'PDE Solver GUI', 'NumberTitle', 'off');

    % Create input fields and labels
    L_text = uicontrol('Style', 'text', 'String', 'Length (L):', 'Position', [50, 320, 100, 20]);
    L_edit = uicontrol('Style', 'edit', 'Position', [150, 320, 100, 20]);

    T_text = uicontrol('Style', 'text', 'String', 'Time (T):', 'Position', [50, 290, 100, 20]);
    T_edit = uicontrol('Style', 'edit', 'Position', [150, 290, 100, 20]);

    Nx_text = uicontrol('Style', 'text', 'String', 'Number of x steps (Nx):', 'Position', [50, 260, 150, 20]);
    Nx_edit = uicontrol('Style', 'edit', 'Position', [200, 260, 100, 20]);

    Nt_text = uicontrol('Style', 'text', 'String', 'Number of time steps (Nt):', 'Position', [50, 230, 150, 20]);
    Nt_edit = uicontrol('Style', 'edit', 'Position', [200, 230, 100, 20]);

    c_text = uicontrol('Style', 'text', 'String', 'Wave speed (c):', 'Position', [50, 200, 100, 20]);
    c_edit = uicontrol('Style', 'edit', 'Position', [150, 200, 100, 20]);

    u0_text = uicontrol('Style', 'text', 'String', 'Initial condition (u0):', 'Position', [50, 170, 150, 20]);
    u0_edit = uicontrol('Style', 'edit', 'Position', [200, 170, 300, 20]);

    du0_dt_text = uicontrol('Style', 'text', 'String', 'Initial derivative (du0_dt):', 'Position', [50, 140, 150, 20]);
    du0_dt_edit = uicontrol('Style', 'edit', 'Position', [200, 140, 300, 20]);

    % Create buttons
    solve_explicit_btn = uicontrol('Style', 'pushbutton', 'String', 'Solve Explicit Method', 'Position', [50, 100, 200, 30], 'Callback', @solveExplicit);
    solve_hyperbolic_btn = uicontrol('Style', 'pushbutton', 'String', 'Solve Hyperbolic PDE', 'Position', [300, 100, 200, 30], 'Callback', @solveHyperbolic);
    sir_model_btn = uicontrol('Style', 'pushbutton', 'String', 'SIR Model 3D', 'Position', [50, 50, 200, 30], 'Callback', @SIRModel3D);
solve_poisson_btn = uicontrol('Style', 'pushbutton', 'String', 'Solve Poisson Equation', 'Position', [300, 50, 200, 30], 'Callback', @solvePoisson);
    % Callback functions
    function solveExplicit(~, ~)
        L = str2double(get(L_edit, 'String'));
        T = str2double(get(T_edit, 'String'));
        Nx = str2double(get(Nx_edit, 'String'));
        Nt = str2double(get(Nt_edit, 'String'));
        lambda = 0.5; % Example value for lambda, you may adjust it as needed
        u = explicitMethod(L, T, Nx, Nt, lambda);
        plotSolution(u, L, T, Nx, Nt);
    end
  function solvePoisson(~, ~)
        % Run the Poisson equation solver script
      Poisson_matrix();

    end
    function solveHyperbolic(~, ~)
        L = str2double(get(L_edit, 'String'));
        T = str2double(get(T_edit, 'String'));
        Nx = str2double(get(Nx_edit, 'String'));
        Nt = str2double(get(Nt_edit, 'String'));
        c = str2double(get(c_edit, 'String'));
        u0_str = get(u0_edit, 'String');
        du0_dt_str = get(du0_dt_edit, 'String');
        
        % Convert input strings to function handles
        u0 = str2func(['@(x)' u0_str]);
        du0_dt = str2func(['@(x)' du0_dt_str]);
        
        u = solveHyperbolicPDE(L, T, Nx, Nt, c, u0, du0_dt);
        plotSolution(u, L, T, Nx, Nt);
    end

    function SIRModel3D(~, ~)
        SIR_model_3D();
    end

    % Plot solution
    function plotSolution(u, L, T, Nx, Nt)
        figure;
        [X, T] = meshgrid(linspace(0, L, Nx), linspace(0, T, Nt + 1));
        surf(X, T, u');
        xlabel('x');
        ylabel('t');
        zlabel('u');
        title('Solution');
    end

    % SIR model function
    function SIR_model_3D()
        % Parameters
        beta = 0.3;     % Infection rate
        gamma = 0.1;    % Recovery rate

        % Initial conditions
        S0 = 0.9;       % Initial proportion of susceptible individuals
        I0 = 0.1;       % Initial proportion of infectious individuals
        R0 = 0;         % Initial proportion of recovered individuals

        y0 = [S0; I0; R0];  % Initial condition vector

        % Time vector
        tspan = [0 100];    % Time span for simulation

        % Solve ODEs
        [t, y] = ode45(@sir_ode, tspan, y0);

        % Plot 3D surface
        figure;
        plot3(y(:,1), y(:,2), y(:,3), 'LineWidth', 2);
        title('SIR Model in 3D');
        xlabel('Susceptible');
        ylabel('Infectious');
        zlabel('Recovered');
        grid on;

        % Define ODEs
        function dydt = sir_ode(t, y)
            S = y(1);   % Susceptible
            I = y(2);   % Infectious
            R = y(3);   % Recovered

            dSdt = -beta * S * I;       % Change in susceptible
            dIdt = beta * S * I - gamma * I;   % Change in infectious
            dRdt = gamma * I;           % Change in recovered

            dydt = [dSdt; dIdt; dRdt];  % Output the changes
        end
    end
end