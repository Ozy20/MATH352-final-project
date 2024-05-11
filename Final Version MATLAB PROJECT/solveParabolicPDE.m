function u = solveParabolicPDE(L, T, Nx, Nt, alpha, u_left, u_right, u0)

dx = L / (Nx - 1);
dt = T / Nt;

x = linspace(0, L, Nx);

A = constructParabolicMatrix(Nx, alpha, dt, dx);

u = zeros(Nx, Nt + 1);
u(:, 1) = u0(x);

for t = 1:Nt
    u(:, t + 1) = A * u(:, t);
    u(1, t + 1) = u_left;
    u(end, t + 1) = u_right;
end

[X, T] = meshgrid(x, linspace(0, T, Nt + 1));
surf(X, T, u');
xlabel('x');
ylabel('t');
zlabel('u');
title('Solution of Parabolic PDE (Heat Equation)');

end

function A = constructParabolicMatrix(Nx, alpha, dt, dx)
A = zeros(Nx, Nx);
A(1, 1) = 1;
A(end, end) = 1;
for i = 2:Nx-1
    A(i, i-1) = -alpha * dt / (dx^2);
    A(i, i) = 1 + 2 * alpha * dt / (dx^2);
    A(i, i+1) = -alpha * dt / (dx^2);
end
end
s
