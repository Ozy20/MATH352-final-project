function u = solveEllipticPDE(L, Nx, Ny, f, u_left, u_right, u_top, u_bottom)

dx = L / (Nx - 1);
dy = L / (Ny - 1);

x = linspace(0, L, Nx);
y = linspace(0, L, Ny);

A = constructEllipticMatrix(Nx, Ny, dx, dy);

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

[X, Y] = meshgrid(x, y);
surf(X, Y, u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Solution of Elliptic PDE (Poisson''s Equation)');

end

function A = constructEllipticMatrix(Nx, Ny, dx, dy)

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

