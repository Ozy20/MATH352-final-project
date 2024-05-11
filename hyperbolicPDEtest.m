L = 1;
T = 0.5;
Nx = 50;
Nt = 100;
c = 1;

u_exact = @(x, t) 0.5 * (sin(pi * (x - c * t)) + sin(pi * (x + c * t)));

u0 = @(x) 0.5 * (sin(pi * x));

du0_dt = @(x) 0;

u = solveHyperbolicPDE(L, T, Nx, Nt, c, u0, du0_dt);

x = linspace(0, L, Nx);
t = linspace(0, T, Nt + 1);
[X, T] = meshgrid(x, t);

u_exact_values = u_exact(X, T);
u_error = abs(u - u_exact_values');

figure;
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


