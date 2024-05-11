classdef app1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        u_rightTextArea_2       matlab.ui.control.TextArea
        u_rightTextArea_2Label  matlab.ui.control.Label
        u_leftTextArea_2        matlab.ui.control.TextArea
        u_leftTextArea_2Label   matlab.ui.control.Label
        alphaTextArea           matlab.ui.control.TextArea
        alphaTextAreaLabel      matlab.ui.control.Label
        LTextArea_3             matlab.ui.control.TextArea
        LTextArea_3Label        matlab.ui.control.Label
        NtTextArea_2            matlab.ui.control.TextArea
        NtTextArea_2Label       matlab.ui.control.Label
        NxTextArea_3            matlab.ui.control.TextArea
        NxTextArea_3Label       matlab.ui.control.Label
        TTextArea_2             matlab.ui.control.TextArea
        TTextArea_2Label        matlab.ui.control.Label
        u_bottomTextArea        matlab.ui.control.TextArea
        u_bottomTextAreaLabel   matlab.ui.control.Label
        u_topTextArea           matlab.ui.control.TextArea
        u_topTextAreaLabel      matlab.ui.control.Label
        u_rightTextArea         matlab.ui.control.TextArea
        u_rightTextAreaLabel    matlab.ui.control.Label
        u_leftTextArea          matlab.ui.control.TextArea
        u_leftTextAreaLabel     matlab.ui.control.Label
        NyTextArea              matlab.ui.control.TextArea
        NyTextAreaLabel         matlab.ui.control.Label
        NxTextArea_2            matlab.ui.control.TextArea
        NxTextArea_2Label       matlab.ui.control.Label
        LTextArea_2             matlab.ui.control.TextArea
        LTextArea_2Label        matlab.ui.control.Label
        cTextArea               matlab.ui.control.TextArea
        cTextAreaLabel          matlab.ui.control.Label
        NtTextArea              matlab.ui.control.TextArea
        NtTextAreaLabel         matlab.ui.control.Label
        NxTextArea              matlab.ui.control.TextArea
        NxTextAreaLabel         matlab.ui.control.Label
        TTextArea               matlab.ui.control.TextArea
        TTextAreaLabel          matlab.ui.control.Label
        LTextArea               matlab.ui.control.TextArea
        LTextAreaLabel          matlab.ui.control.Label
        DropDown                matlab.ui.control.DropDown
        DropDownLabel           matlab.ui.control.Label
        PlotButton              matlab.ui.control.Button
    end

    
    
    methods (Access = private)
       function u = solveEllipticPDE(app,L, Nx, Ny, f, u_left, u_right, u_top, u_bottom)

dx = L / (Nx - 1);
dy = L / (Ny - 1);

x = linspace(0, L, Nx);
y = linspace(0, L, Ny);

A = constructEllipticMatrix(app,Nx, Ny, dx, dy);

F = zeros(Nx * Ny, 1);
for i = 1:Nx
    for j = 1:Ny
        idx = (i - 1) * Ny + j;
        F(idx) = f(x(i), y(j));
    end
end

u = A \ F;

u = reshape(u, Ny, Nx)';

u(:, 1) = u_left;
u(:, end) = u_right;
u(1, :) = u_bottom;
u(end, :) = u_top;
figure('Name', '2')
[X, Y] = meshgrid(x, y);
surf(X, Y, u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Solution of Elliptic PDE (Poisson''s Equation)');

end

function A = constructEllipticMatrix(app,Nx, Ny, dx, dy)

N = Nx * Ny;
A = zeros(N, N);

for i = 1:N
    A(i, i) = -2 * (1 / dx^2 + 1 / dy^2);
    if i - 1 > 0 && mod(i - 1, Ny) ~= 0
        A(i, i - 1) = 1 / dx^2;
    end
    if i + 1 <= N && mod(i, Ny) ~= 0
        A(i, i + 1) = 1 / dx^2;
    end
    if i - Ny > 0
        A(i, i - Ny) = 1 / dy^2;
    end
    if i + Ny <= N
        A(i, i + Ny) = 1 / dy^2;
    end
end
end 
        
        
        function u = solveParabolicPDE(app,L, T, Nx, Nt, alpha, u_left, u_right, u0)

dx = L / (Nx - 1);
dt = T / Nt;

x = linspace(0, L, Nx);

A = constructParabolicMatrix(app,Nx, alpha, dt, dx);

u = zeros(Nx, Nt + 1);
u(:, 1) = u0(x);

for t = 1:Nt
    u(:, t + 1) = A * u(:, t);
    u(1, t + 1) = u_left;
    u(end, t + 1) = u_right;
end
figure('Name', '2')
[X, T] = meshgrid(x, linspace(0, T, Nt + 1));
surf(X, T, u');
xlabel('x');
ylabel('t');
zlabel('u');
title('Solution of Parabolic PDE (Heat Equation)');

end

function A = constructParabolicMatrix(app,Nx, alpha, dt, dx)
A = zeros(Nx, Nx);
A(1, 1) = 1;
A(end, end) = 1;
for i = 2:Nx-1
    A(i, i-1) = -alpha * dt / (dx^2);
    A(i, i) = 1 + 2 * alpha * dt / (dx^2);
    A(i, i+1) = -alpha * dt / (dx^2);
end
end

      function [t, y] = solveODESystem(app,f, tspan, y0, N)
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
      end;
        
        function u = solveHyperbolicPDE(app,L, T, Nx, Nt, c, u0, du0_dt)

dx = L / (Nx - 1);
dt = T / Nt;

x = linspace(0, L, Nx);

A = constructHyperbolicMatrix(app,Nx, c, dt, dx);

u = zeros(Nx, Nt + 1);
u(:, 1) = u0(x);

u(:, 2) = u(:, 1) + dt * du0_dt(x);

for t = 2:Nt
    u(:, t + 1) = A * u(:, t) - u(:, t - 1);
end
figure('Name', '2')
[X, T] = meshgrid(x, linspace(0, T, Nt + 1));
surf(X, T, u');
xlabel('x');
ylabel('t');
zlabel('u');
title('Solution of Hyperbolic PDE (Wave Equation)');

end

function A = constructHyperbolicMatrix(app,Nx, c, dt, dx)

A = eye(Nx);
A = A + (c * dt / dx) * circshift(eye(Nx), -1, 2);
A = A - (c * dt / dx) * circshift(eye(Nx), 1, 2);
end
    end;


  
        
    


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
           if(app.DropDown.Value=="Poisson")
            L = str2double(app.LTextArea_2.Value);
Nx = str2double(app.NxTextArea_2.Value);
Ny = str2double(app.NyTextArea.Value);
u_left = str2double(app.u_leftTextArea.Value);
u_right = str2double(app.u_rightTextArea.Value);
u_top = str2double(app.u_topTextArea.Value);
u_bottom = str2double(app.u_bottomTextArea.Value);

u_exact = @(x, y) sin(pi * x) * sin(pi * y);

f = @(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y);

u = solveEllipticPDE(app,L, Nx, Ny, f, u_left, u_right, u_top, u_bottom);

[x, y] = meshgrid(linspace(0, L, Nx), linspace(0, L, Ny));
u_error = abs(u - u_exact(x, y));

figure('Name', '1')
subplot(1, 2, 1);
surf(x, y, u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Numerical Solution');

subplot(1, 2, 2);
surf(x, y, u_error);
xlabel('x');
ylabel('y');
zlabel('Error');
title('Error');
           else if(app.DropDown.Value=="Heat")
                  L = str2double(app.LTextArea_3.Value);
T = str2double(app.TTextArea_2.Value);
Nx = str2double(app.NxTextArea_3.Value);
Nt = str2double(app.NtTextArea_2.Value); 
alpha = str2double(app.alphaTextArea.Value);
u_left = str2double(app.u_leftTextArea_2.Value);
u_right = str2double(app.u_rightTextArea_2.Value);

u_exact = @(x, t) exp(-pi^2 * alpha * t') * sin(pi * x);

u0 = @(x) sin(pi * x);

u = solveParabolicPDE(app,L, T, Nx, Nt, alpha, u_left, u_right, u0);

x = linspace(0, L, Nx);
t = linspace(0, T, Nt + 1);
u_exact_values = u_exact(x, t);

u_exact_values = u_exact_values';

u_error = abs(u - u_exact_values);

figure('Name', '1')
subplot(1, 2, 1);
surf(t, x, u);
xlabel('t');
ylabel('x');
zlabel('u');
title('Numerical Solution');

subplot(1, 2, 2);
surf(t, x, u_error);
xlabel('t');
ylabel('x');
zlabel('Error');
title('Error');
           else if(app.DropDown.Value=="systemofOdes")
                   
                   % Define the system of ODEs
% Example: dy1/dt = -2*y1 + y2, dy2/dt = y1 - 2*y2, dy3/dt = -y1 + y3
f = @(t, y) [-2*y(1) + y(2); y(1) - 2*y(2); -y(1) + y(3)];

% Define the time span and initial conditions
tspan = [0, 2];
y0 = [1; 0; 2]; % Initial conditions for y1, y2, and y3

% Number of time steps
N = 100;

% Solve the system of ODEs using RK4 method
[t, y] = solveODESystem(app,f, tspan, y0, N);

% Plot the solution
figure('Name', '1')
subplot(3, 1, 1);
plot(t, y(1, :), '-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('y_1(t)');
title('Solution of dy_1/dt = -2*y_1 + y_2');

subplot(3, 1, 2);
plot(t, y(2, :), '-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('y_2(t)');
title('Solution of dy_2/dt = y_1 - 2*y_2');
subplot(3, 1, 3);
plot(t, y(3, :), '-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('y_3(t)');
title('Solution of dy_3/dt = -y_1 + y_3');
           else
               L = str2double(app.LTextArea.Value);
T = str2double(app.TTextArea.Value);
Nx = str2double(app.NxTextArea.Value);
Nt = str2double(app.NtTextArea.Value);
c = str2double(app.cTextArea.Value);

u_exact = @(x, t) 0.5 * (sin(pi * (x - c * t)) + sin(pi * (x + c * t)));

u0 = @(x) 0.5 * (sin(pi * x));

du0_dt = @(x) 0;

u = solveHyperbolicPDE(app,L, T, Nx, Nt, c, u0, du0_dt);

x = linspace(0, L, Nx);
t = linspace(0, T, Nt + 1);
[X, T] = meshgrid(x, t);

u_exact_values = u_exact(X, T);
u_error = abs(u - u_exact_values');

figure('Name', '1')
subplot(1, 2, 1);
surf(X, T, u');
xlabel('x');
ylabel('t');
zlabel('u');
title('Numerical Solution');

subplot(1, 2, 2);
surf(X, T, u_error');
xlabel('x');
ylabel('t');
zlabel('Error');
title('Error');
           end;       
           end;
           end;
        
        end

        % Value changed function: DropDown
        function DropDownValueChanged(app, event)
            if(app.DropDown.Value=="Wave")




                app.TTextArea.Visible='on';
                  app.NxTextArea.Visible='on';
                  app.NtTextArea.Visible='on';
                  app.LTextArea.Visible='on';
                  app.cTextArea.Visible='on';

                  app.TTextAreaLabel.Visible='on';
                  app.NxTextAreaLabel.Visible='on';
                  app.NtTextAreaLabel.Visible='on';
                  app.LTextAreaLabel.Visible='on';
                  app.cTextAreaLabel.Visible='on';







                app.LTextArea_2.Visible='off';
                app.NxTextArea_2.Visible='off';
               app.NyTextArea.Visible='off';
                app.u_leftTextArea.Visible='off';
                app.u_rightTextArea.Visible='off';
              app.u_bottomTextArea.Visible='off';
              app.u_topTextArea.Visible='off';
              app.LTextArea_2Label.Visible='off';
                app.NxTextArea_2Label.Visible='off';
               app.NyTextAreaLabel.Visible='off';
                app.u_leftTextAreaLabel.Visible='off';
                app.u_rightTextAreaLabel.Visible='off';
              app.u_bottomTextAreaLabel.Visible='off';
              app.u_topTextAreaLabel.Visible='off';

              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

              app.LTextArea_3.Visible='off';
                app.NxTextArea_3.Visible='off';
                app.NtTextArea_2.Visible='off';
                app.TTextArea_2.Visible='off';
                app.u_leftTextArea_2.Visible='off';
                app.u_rightTextArea_2.Visible='off';
                app.alphaTextArea.Visible='off';

                app.LTextArea_3Label.Visible='off';
                app.NxTextArea_3Label.Visible='off';
                app.NtTextArea_2Label.Visible='off';
                app.TTextArea_2Label.Visible='off';
                app.u_leftTextArea_2Label.Visible='off';
                app.u_rightTextArea_2Label.Visible='off';
                app.alphaTextAreaLabel.Visible='off';


            else if(app.DropDown.Value=="Heat")


                    app.LTextArea_3.Visible='on';
                app.NxTextArea_3.Visible='on';
                app.NtTextArea_2.Visible='on';
                app.TTextArea_2.Visible='on';
                app.u_leftTextArea_2.Visible='on';
                app.u_rightTextArea_2.Visible='on';
                app.alphaTextArea.Visible='on';

                app.LTextArea_3Label.Visible='on';
                app.NxTextArea_3Label.Visible='on';
                app.NtTextArea_2Label.Visible='on';
                app.TTextArea_2Label.Visible='on';
                app.u_leftTextArea_2Label.Visible='on';
                app.u_rightTextArea_2Label.Visible='on';
                app.alphaTextAreaLabel.Visible='on';




                 app.LTextArea_2.Visible='off';
                app.NxTextArea_2.Visible='off';
               app.NyTextArea.Visible='off';
                app.u_leftTextArea.Visible='off';
                app.u_rightTextArea.Visible='off';
              app.u_bottomTextArea.Visible='off';
              app.u_topTextArea.Visible='off';
              app.LTextArea_2Label.Visible='off';
                app.NxTextArea_2Label.Visible='off';
               app.NyTextAreaLabel.Visible='off';
                app.u_leftTextAreaLabel.Visible='off';
                app.u_rightTextAreaLabel.Visible='off';
              app.u_bottomTextAreaLabel.Visible='off';
              app.u_topTextAreaLabel.Visible='off';

              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  app.TTextArea.Visible='off';
                  app.NxTextArea.Visible='off';
                  app.NtTextArea.Visible='off';
                  app.LTextArea.Visible='off';
                  app.cTextArea.Visible='off';

                  app.TTextAreaLabel.Visible='off';
                  app.NxTextAreaLabel.Visible='off';
                  app.NtTextAreaLabel.Visible='off';
                  app.LTextAreaLabel.Visible='off';
                  app.cTextAreaLabel.Visible='off';

                  else if(app.DropDown.Value=="Poisson")

                          app.LTextArea_2.Visible='on';
                app.NxTextArea_2.Visible='on';
               app.NyTextArea.Visible='on';
                app.u_leftTextArea.Visible='on';
                app.u_rightTextArea.Visible='on';
              app.u_bottomTextArea.Visible='on';
              app.u_topTextArea.Visible='on';
              app.LTextArea_2Label.Visible='on';
                app.NxTextArea_2Label.Visible='on';
               app.NyTextAreaLabel.Visible='on';
                app.u_leftTextAreaLabel.Visible='on';
                app.u_rightTextAreaLabel.Visible='on';
              app.u_bottomTextAreaLabel.Visible='on';
              app.u_topTextAreaLabel.Visible='on';




                          app.TTextArea.Visible='off';
                  app.NxTextArea.Visible='off';
                  app.NtTextArea.Visible='off';
                  app.LTextArea.Visible='off';
                  app.cTextArea.Visible='off';

                  app.TTextAreaLabel.Visible='off';
                  app.NxTextAreaLabel.Visible='off';
                  app.NtTextAreaLabel.Visible='off';
                  app.LTextAreaLabel.Visible='off';
                  app.cTextAreaLabel.Visible='off';

                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                  app.LTextArea_3.Visible='off';
                app.NxTextArea_3.Visible='off';
                app.NtTextArea_2.Visible='off';
                app.TTextArea_2.Visible='off';
                app.u_leftTextArea_2.Visible='off';
                app.u_rightTextArea_2.Visible='off';
                app.alphaTextArea.Visible='off';

                app.LTextArea_3Label.Visible='off';
                app.NxTextArea_3Label.Visible='off';
                app.NtTextArea_2Label.Visible='off';
                app.TTextArea_2Label.Visible='off';
                app.u_leftTextArea_2Label.Visible='off';
                app.u_rightTextArea_2Label.Visible='off';
                app.alphaTextAreaLabel.Visible='off';

                  else

                      app.LTextArea_3.Visible='off';
                app.NxTextArea_3.Visible='off';
                app.NtTextArea_2.Visible='off';
                app.TTextArea_2.Visible='off';
                app.u_leftTextArea_2.Visible='off';
                app.u_rightTextArea_2.Visible='off';
                app.alphaTextArea.Visible='off';

                app.LTextArea_3Label.Visible='off';
                app.NxTextArea_3Label.Visible='off';
                app.NtTextArea_2Label.Visible='off';
                app.TTextArea_2Label.Visible='off';
                app.u_leftTextArea_2Label.Visible='off';
                app.u_rightTextArea_2Label.Visible='off';
                app.alphaTextAreaLabel.Visible='off';

                 app.TTextArea.Visible='off';
                  app.NxTextArea.Visible='off';
                  app.NtTextArea.Visible='off';
                  app.LTextArea.Visible='off';
                  app.cTextArea.Visible='off';

                  app.TTextAreaLabel.Visible='off';
                  app.NxTextAreaLabel.Visible='off';
                  app.NtTextAreaLabel.Visible='off';
                  app.LTextAreaLabel.Visible='off';
                  app.cTextAreaLabel.Visible='off';

                  app.LTextArea_2.Visible='off';
                app.NxTextArea_2.Visible='off';
               app.NyTextArea.Visible='off';
                app.u_leftTextArea.Visible='off';
                app.u_rightTextArea.Visible='off';
              app.u_bottomTextArea.Visible='off';
              app.u_topTextArea.Visible='off';
              app.LTextArea_2Label.Visible='off';
                app.NxTextArea_2Label.Visible='off';
               app.NyTextAreaLabel.Visible='off';
                app.u_leftTextAreaLabel.Visible='off';
                app.u_rightTextAreaLabel.Visible='off';
              app.u_bottomTextAreaLabel.Visible='off';
              app.u_topTextAreaLabel.Visible='off';



                  end
            end
            end





              
            
            

  % Additional customization (optional)
  
        
    

            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.902 0.902 0.902];
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create PlotButton
            app.PlotButton = uibutton(app.UIFigure, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.Position = [271 56 100 22];
            app.PlotButton.Text = 'Plot';

            % Create DropDownLabel
            app.DropDownLabel = uilabel(app.UIFigure);
            app.DropDownLabel.HorizontalAlignment = 'right';
            app.DropDownLabel.FontColor = [0.149 0.149 0.149];
            app.DropDownLabel.Position = [288 155 65 22];
            app.DropDownLabel.Text = 'Drop Down';

            % Create DropDown
            app.DropDown = uidropdown(app.UIFigure);
            app.DropDown.Items = {'Heat', 'Poisson', 'systemofOdes', 'Wave'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            app.DropDown.Position = [271 97 100 22];
            app.DropDown.Value = 'systemofOdes';

            % Create LTextAreaLabel
            app.LTextAreaLabel = uilabel(app.UIFigure);
            app.LTextAreaLabel.HorizontalAlignment = 'right';
            app.LTextAreaLabel.Visible = 'off';
            app.LTextAreaLabel.Position = [4 444 25 22];
            app.LTextAreaLabel.Text = 'L';

            % Create LTextArea
            app.LTextArea = uitextarea(app.UIFigure);
            app.LTextArea.Visible = 'off';
            app.LTextArea.Position = [44 444 131 24];

            % Create TTextAreaLabel
            app.TTextAreaLabel = uilabel(app.UIFigure);
            app.TTextAreaLabel.HorizontalAlignment = 'right';
            app.TTextAreaLabel.Visible = 'off';
            app.TTextAreaLabel.Position = [10 398 25 22];
            app.TTextAreaLabel.Text = 'T';

            % Create TTextArea
            app.TTextArea = uitextarea(app.UIFigure);
            app.TTextArea.Visible = 'off';
            app.TTextArea.Position = [50 398 131 24];

            % Create NxTextAreaLabel
            app.NxTextAreaLabel = uilabel(app.UIFigure);
            app.NxTextAreaLabel.HorizontalAlignment = 'right';
            app.NxTextAreaLabel.Visible = 'off';
            app.NxTextAreaLabel.Position = [9 352 25 22];
            app.NxTextAreaLabel.Text = 'Nx';

            % Create NxTextArea
            app.NxTextArea = uitextarea(app.UIFigure);
            app.NxTextArea.Visible = 'off';
            app.NxTextArea.Position = [49 352 131 24];

            % Create NtTextAreaLabel
            app.NtTextAreaLabel = uilabel(app.UIFigure);
            app.NtTextAreaLabel.HorizontalAlignment = 'right';
            app.NtTextAreaLabel.Visible = 'off';
            app.NtTextAreaLabel.Position = [10 303 25 22];
            app.NtTextAreaLabel.Text = 'Nt';

            % Create NtTextArea
            app.NtTextArea = uitextarea(app.UIFigure);
            app.NtTextArea.Visible = 'off';
            app.NtTextArea.Position = [50 303 131 24];

            % Create cTextAreaLabel
            app.cTextAreaLabel = uilabel(app.UIFigure);
            app.cTextAreaLabel.HorizontalAlignment = 'right';
            app.cTextAreaLabel.Visible = 'off';
            app.cTextAreaLabel.Position = [14 251 25 22];
            app.cTextAreaLabel.Text = 'c';

            % Create cTextArea
            app.cTextArea = uitextarea(app.UIFigure);
            app.cTextArea.Visible = 'off';
            app.cTextArea.Position = [54 251 131 24];

            % Create LTextArea_2Label
            app.LTextArea_2Label = uilabel(app.UIFigure);
            app.LTextArea_2Label.HorizontalAlignment = 'right';
            app.LTextArea_2Label.Visible = 'off';
            app.LTextArea_2Label.Position = [4 444 25 22];
            app.LTextArea_2Label.Text = 'L';

            % Create LTextArea_2
            app.LTextArea_2 = uitextarea(app.UIFigure);
            app.LTextArea_2.Visible = 'off';
            app.LTextArea_2.Position = [44 444 131 24];

            % Create NxTextArea_2Label
            app.NxTextArea_2Label = uilabel(app.UIFigure);
            app.NxTextArea_2Label.HorizontalAlignment = 'right';
            app.NxTextArea_2Label.Visible = 'off';
            app.NxTextArea_2Label.Position = [10 397 25 22];
            app.NxTextArea_2Label.Text = 'Nx';

            % Create NxTextArea_2
            app.NxTextArea_2 = uitextarea(app.UIFigure);
            app.NxTextArea_2.Visible = 'off';
            app.NxTextArea_2.Position = [50 397 131 24];

            % Create NyTextAreaLabel
            app.NyTextAreaLabel = uilabel(app.UIFigure);
            app.NyTextAreaLabel.HorizontalAlignment = 'right';
            app.NyTextAreaLabel.Visible = 'off';
            app.NyTextAreaLabel.Position = [10 350 25 22];
            app.NyTextAreaLabel.Text = 'Ny';

            % Create NyTextArea
            app.NyTextArea = uitextarea(app.UIFigure);
            app.NyTextArea.Visible = 'off';
            app.NyTextArea.Position = [50 350 131 24];

            % Create u_leftTextAreaLabel
            app.u_leftTextAreaLabel = uilabel(app.UIFigure);
            app.u_leftTextAreaLabel.HorizontalAlignment = 'right';
            app.u_leftTextAreaLabel.Visible = 'off';
            app.u_leftTextAreaLabel.Position = [1 302 34 22];
            app.u_leftTextAreaLabel.Text = 'u_left';

            % Create u_leftTextArea
            app.u_leftTextArea = uitextarea(app.UIFigure);
            app.u_leftTextArea.Visible = 'off';
            app.u_leftTextArea.Position = [50 302 131 24];

            % Create u_rightTextAreaLabel
            app.u_rightTextAreaLabel = uilabel(app.UIFigure);
            app.u_rightTextAreaLabel.HorizontalAlignment = 'right';
            app.u_rightTextAreaLabel.Visible = 'off';
            app.u_rightTextAreaLabel.Position = [-3 252 42 22];
            app.u_rightTextAreaLabel.Text = 'u_right';

            % Create u_rightTextArea
            app.u_rightTextArea = uitextarea(app.UIFigure);
            app.u_rightTextArea.Visible = 'off';
            app.u_rightTextArea.Position = [54 252 131 24];

            % Create u_topTextAreaLabel
            app.u_topTextAreaLabel = uilabel(app.UIFigure);
            app.u_topTextAreaLabel.HorizontalAlignment = 'right';
            app.u_topTextAreaLabel.Visible = 'off';
            app.u_topTextAreaLabel.Position = [1 211 35 22];
            app.u_topTextAreaLabel.Text = 'u_top';

            % Create u_topTextArea
            app.u_topTextArea = uitextarea(app.UIFigure);
            app.u_topTextArea.Visible = 'off';
            app.u_topTextArea.Position = [51 211 131 24];

            % Create u_bottomTextAreaLabel
            app.u_bottomTextAreaLabel = uilabel(app.UIFigure);
            app.u_bottomTextAreaLabel.HorizontalAlignment = 'right';
            app.u_bottomTextAreaLabel.Visible = 'off';
            app.u_bottomTextAreaLabel.Position = [4 167 55 22];
            app.u_bottomTextAreaLabel.Text = 'u_bottom';

            % Create u_bottomTextArea
            app.u_bottomTextArea = uitextarea(app.UIFigure);
            app.u_bottomTextArea.Visible = 'off';
            app.u_bottomTextArea.Position = [74 167 131 24];

            % Create TTextArea_2Label
            app.TTextArea_2Label = uilabel(app.UIFigure);
            app.TTextArea_2Label.HorizontalAlignment = 'right';
            app.TTextArea_2Label.Visible = 'off';
            app.TTextArea_2Label.Position = [10 398 25 22];
            app.TTextArea_2Label.Text = 'T';

            % Create TTextArea_2
            app.TTextArea_2 = uitextarea(app.UIFigure);
            app.TTextArea_2.Visible = 'off';
            app.TTextArea_2.Position = [50 398 131 24];

            % Create NxTextArea_3Label
            app.NxTextArea_3Label = uilabel(app.UIFigure);
            app.NxTextArea_3Label.HorizontalAlignment = 'right';
            app.NxTextArea_3Label.Visible = 'off';
            app.NxTextArea_3Label.Position = [10 351 25 22];
            app.NxTextArea_3Label.Text = 'Nx';

            % Create NxTextArea_3
            app.NxTextArea_3 = uitextarea(app.UIFigure);
            app.NxTextArea_3.Visible = 'off';
            app.NxTextArea_3.Position = [50 351 131 24];

            % Create NtTextArea_2Label
            app.NtTextArea_2Label = uilabel(app.UIFigure);
            app.NtTextArea_2Label.HorizontalAlignment = 'right';
            app.NtTextArea_2Label.Visible = 'off';
            app.NtTextArea_2Label.Position = [9 303 25 22];
            app.NtTextArea_2Label.Text = 'Nt';

            % Create NtTextArea_2
            app.NtTextArea_2 = uitextarea(app.UIFigure);
            app.NtTextArea_2.Visible = 'off';
            app.NtTextArea_2.Position = [49 303 131 24];

            % Create LTextArea_3Label
            app.LTextArea_3Label = uilabel(app.UIFigure);
            app.LTextArea_3Label.HorizontalAlignment = 'right';
            app.LTextArea_3Label.Visible = 'off';
            app.LTextArea_3Label.Position = [4 444 25 22];
            app.LTextArea_3Label.Text = 'L';

            % Create LTextArea_3
            app.LTextArea_3 = uitextarea(app.UIFigure);
            app.LTextArea_3.Visible = 'off';
            app.LTextArea_3.Position = [44 444 131 24];

            % Create alphaTextAreaLabel
            app.alphaTextAreaLabel = uilabel(app.UIFigure);
            app.alphaTextAreaLabel.HorizontalAlignment = 'right';
            app.alphaTextAreaLabel.Visible = 'off';
            app.alphaTextAreaLabel.Position = [0 254 34 22];
            app.alphaTextAreaLabel.Text = 'alpha';

            % Create alphaTextArea
            app.alphaTextArea = uitextarea(app.UIFigure);
            app.alphaTextArea.Visible = 'off';
            app.alphaTextArea.Position = [49 254 131 24];

            % Create u_leftTextArea_2Label
            app.u_leftTextArea_2Label = uilabel(app.UIFigure);
            app.u_leftTextArea_2Label.HorizontalAlignment = 'right';
            app.u_leftTextArea_2Label.Visible = 'off';
            app.u_leftTextArea_2Label.Position = [5 212 34 22];
            app.u_leftTextArea_2Label.Text = 'u_left';

            % Create u_leftTextArea_2
            app.u_leftTextArea_2 = uitextarea(app.UIFigure);
            app.u_leftTextArea_2.Visible = 'off';
            app.u_leftTextArea_2.Position = [54 212 131 24];

            % Create u_rightTextArea_2Label
            app.u_rightTextArea_2Label = uilabel(app.UIFigure);
            app.u_rightTextArea_2Label.HorizontalAlignment = 'right';
            app.u_rightTextArea_2Label.Visible = 'off';
            app.u_rightTextArea_2Label.Position = [5 168 42 22];
            app.u_rightTextArea_2Label.Text = 'u_right';

            % Create u_rightTextArea_2
            app.u_rightTextArea_2 = uitextarea(app.UIFigure);
            app.u_rightTextArea_2.Visible = 'off';
            app.u_rightTextArea_2.Position = [62 167 131 24];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end