function u = solveHyperbolicPDE(L, T, Nx, Nt, c, u0, du0_dt)

dx = L / (Nx - 1);
dt = T / Nt;

x = linspace(0, L, Nx);

A = constructHyperbolicMatrix(Nx, c, dt, dx);

u = zeros(Nx, Nt + 1);
u(:, 1) = u0(x);

u(:, 2) = u(:, 1) + dt * du0_dt(x);

for t = 2:Nt
    u(:, t + 1) = A * u(:, t) - u(:, t - 1);
end

[X, T] = meshgrid(x, linspace(0, T, Nt + 1));
surf(X, T, u');
xlabel('x');
ylabel('t');
zlabel('u');
title('Solution of Hyperbolic PDE (Wave Equation)');

end

function A = constructHyperbolicMatrix(Nx, c, dt, dx)

A = eye(Nx);
A = A + (c * dt / dx) * circshift(eye(Nx), -1, 2);
A = A - (c * dt / dx) * circshift(eye(Nx), 1, 2);
end

