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
