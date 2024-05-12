function u = explicitMethod(L, T, Nx, Nt, lambda)
    dx = L / (Nx - 1);
    x = 0:dx:L;

    dt = T / Nt;
    t = 0:dt:T;

    u = zeros(Nx, Nt+1);

    % Set the initial condition
    for i = 1:Nx
        u(i,1) = sin(pi*x(i));
    end

    for n = 1:Nt
        % Left boundary (i = 1)
        u(1,n+1) = u(1,n) + lambda * (u(2,n) - u(1,n)) * dt;

        % Interior points
        for i = 2:Nx-1
            u(i,n+1) = u(i,n) + lambda * (u(i+1,n) - 2*u(i,n) + u(i-1,n)) * dt;
        end

        % Right boundary (i = Nx)
        u(Nx,n+1) = u(Nx,n) + lambda * (u(Nx-1,n) - u(Nx,n)) * dt;
    end
end


