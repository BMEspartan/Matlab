%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lid_driven_cavity_navier_stokes.m
%
% Educational 2D incompressible Navier–Stokes solver for a lid-driven
% cavity using a projection method and finite differences.
%
%   ∂u/∂t + (u·∇)u = -∇p + (1/Re) ∇²u
%   ∇·u = 0
%
% Domain: [0,1] x [0,1]
% BCs:
%   - No-slip on left, right, bottom walls (u = v = 0)
%   - Moving lid on top wall (u = U_lid, v = 0)
%
% This script:
%   - Discretizes the equations on a uniform Cartesian grid
%   - Uses explicit time stepping for convection & diffusion
%   - Solves a pressure Poisson equation each step to enforce ∇·u = 0
%   - Visualizes the final velocity field and streamlines
%
% Author: Brian Lindskog
% GitHub: https://github.com/BMEspartan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% ------------------------- Problem Parameters ------------------------ %%
Re      = 100;         % Reynolds number (dimensionless)
U_lid   = 1.0;         % Lid velocity (non-dimensional)
rho     = 1.0;         % Density (non-dimensional)
Lx      = 1.0;         % Domain length in x
Ly      = 1.0;         % Domain length in y

Nx      = 41;          % Number of grid points in x
Ny      = 41;          % Number of grid points in y

dx      = Lx/(Nx-1);
dy      = Ly/(Ny-1);

x       = linspace(0,Lx,Nx);
y       = linspace(0,Ly,Ny);
[X,Y]   = meshgrid(x,y);   % X(i,j): x, Y(i,j): y

% Time-stepping parameters
dt      = 0.001;       % Time step (adjust for stability / speed)
nsteps  = 2000;        % Number of time steps (t_final ~ nsteps*dt)
nit     = 50;          % Pressure Poisson iterations per step

% Visualization frequency
plotEvery = 200;       % Plot every "plotEvery" time steps

%% -------------------------- Field Variables -------------------------- %%
% u(i,j) and v(i,j) are velocities at grid nodes (collocated grid)
u = zeros(Ny,Nx);      % u-velocity
v = zeros(Ny,Nx);      % v-velocity
p = zeros(Ny,Nx);      % pressure
b = zeros(Ny,Nx);      % RHS for pressure Poisson equation

%% ------------------------ Helper: Boundary Conditions ---------------- %%
applyVelocityBC = @() ...
    ( ...
        % No-slip on left and right walls
        (u(:,1)  = 0)       & (v(:,1)  = 0) & ...
        (u(:,end)= 0)       & (v(:,end)= 0) & ...
        % No-slip on bottom wall
        (u(1,:)  = 0)       & (v(1,:)  = 0) & ...
        % Moving lid on top wall
        (u(end,:)= U_lid)   & (v(end,:)= 0) ...
    );

applyPressureBC = @() ...
    ( ...
        % Zero-gradient (Neumann) on all walls except one reference
        (p(:,end) = p(:,end-1)) & ...  % top
        (p(:,1)   = p(:,2))     & ...  % bottom
        (p(1,:)   = p(2,:))     & ...  % left
        (p(end,:) = 0)                % right wall as reference (p=0)
    );

%% -------------------------- Time Integration -------------------------- %%
fprintf('Starting lid-driven cavity simulation...\n');
fprintf('Grid: %d x %d, Re = %.1f, dt = %.4f, steps = %d\n',...
    Nx, Ny, Re, dt, nsteps);

for n = 1:nsteps
    % Store previous time step
    un = u;
    vn = v;
    pn = p;

    %% ---- Build RHS of pressure Poisson equation: b ----
    % b ≈ ρ [ (1/Δt)(∂u/∂x + ∂v/∂y)
    %       - (∂u/∂x)^2 - 2(∂u/∂y ∂v/∂x) - (∂v/∂y)^2 ]
    dudx = (un(:,3:end)   - un(:,1:end-2)) / (2*dx);  % Ny x (Nx-2)
    dvdy = (vn(3:end,:)   - vn(1:end-2,:)) / (2*dy);  % (Ny-2) x Nx

    dudy = (un(3:end,:)   - un(1:end-2,:)) / (2*dy);  % (Ny-2) x Nx
    dvdx = (vn(:,3:end)   - vn(:,1:end-2)) / (2*dx);  % Ny x (Nx-2)

    % We’ll index interior points (2:Ny-1, 2:Nx-1) carefully:
    b(2:Ny-1,2:Nx-1) = rho * ( ...
        (1/dt) * ( ...
            (un(2:Ny-1,3:Nx) - un(2:Ny-1,1:Nx-2)) / (2*dx) + ...
            (vn(3:Ny,2:Nx-1) - vn(1:Ny-2,2:Nx-1)) / (2*dy) ) ...
        - ...
          ((un(2:Ny-1,3:Nx) - un(2:Ny-1,1:Nx-2)) / (2*dx)).^2 ...
        - 2 * ( ...
            ( (un(3:Ny,2:Nx-1) - un(1:Ny-2,2:Nx-1)) / (2*dy) ) .* ...
            ( (vn(2:Ny-1,3:Nx) - vn(2:Ny-1,1:Nx-2)) / (2*dx) ) ) ...
        - ...
          ((vn(3:Ny,2:Nx-1) - vn(1:Ny-2,2:Nx-1)) / (2*dy)).^2 ...
    );

    %% ---- Pressure Poisson solver (Gauss–Seidel iterations) ----
    p = pn;  % start from previous pressure

    for it = 1:nit
        p_old = p;

        p(2:Ny-1,2:Nx-1) = ( ...
            (p_old(2:Ny-1,3:Nx) + p_old(2:Ny-1,1:Nx-2)) * dy^2 + ...
            (p_old(3:Ny,2:Nx-1) + p_old(1:Ny-2,2:Nx-1)) * dx^2 - ...
            b(2:Ny-1,2:Nx-1) * dx^2 * dy^2 ...
        ) / (2 * (dx^2 + dy^2));

        applyPressureBC();
    end

    %% ---- Velocity update: explicit convection-diffusion + pressure ----
    % u-momentum
    u(2:Ny-1,2:Nx-1) = un(2:Ny-1,2:Nx-1) ...
        - dt * ( ...
            un(2:Ny-1,2:Nx-1) .* ...
            (un(2:Ny-1,3:Nx) - un(2:Ny-1,1:Nx-2)) / (2*dx) + ...
            vn(2:Ny-1,2:Nx-1) .* ...
            (un(3:Ny,2:Nx-1) - un(1:Ny-2,2:Nx-1)) / (2*dy) ...
        ) ...
        - dt / (2*rho*dx) * (p(2:Ny-1,3:Nx) - p(2:Ny-1,1:Nx-2)) ...
        + dt * (1/Re) * ( ...
            (un(2:Ny-1,3:Nx) - 2*un(2:Ny-1,2:Nx-1) + ...
             un(2:Ny-1,1:Nx-2)) / dx^2 + ...
            (un(3:Ny,2:Nx-1) - 2*un(2:Ny-1,2:Nx-1) + ...
             un(1:Ny-2,2:Nx-1)) / dy^2 ...
        );

    % v-momentum
    v(2:Ny-1,2:Nx-1) = vn(2:Ny-1,2:Nx-1) ...
        - dt * ( ...
            un(2:Ny-1,2:Nx-1) .* ...
            (vn(2:Ny-1,3:Nx) - vn(2:Ny-1,1:Nx-2)) / (2*dx) + ...
            vn(2:Ny-1,2:Nx-1) .* ...
            (vn(3:Ny,2:Nx-1) - vn(1:Ny-2,2:Nx-1)) / (2*dy) ...
        ) ...
        - dt / (2*rho*dy) * (p(3:Ny,2:Nx-1) - p(1:Ny-2,2:Nx-1)) ...
        + dt * (1/Re) * ( ...
            (vn(2:Ny-1,3:Nx) - 2*vn(2:Ny-1,2:Nx-1) + ...
             vn(2:Ny-1,1:Nx-2)) / dx^2 + ...
            (vn(3:Ny,2:Nx-1) - 2*vn(2:Ny-1,2:Nx-1) + ...
             vn(1:Ny-2,2:Nx-1)) / dy^2 ...
        );

    % Re-apply velocity boundary conditions after update
    applyVelocityBC();

    %% ------------------------- Visualization ------------------------ %%
    if mod(n, plotEvery) == 0 || n == nsteps
        speed = sqrt(u.^2 + v.^2);

        figure(1); clf;
        contourf(X, Y, speed, 20, 'LineStyle', 'none');
        hold on;
        quiver(X, Y, u, v, 2, 'k');
        hold off;
        colorbar;
        xlabel('x'); ylabel('y');
        title(sprintf('Lid-Driven Cavity (Re = %.1f), t = %.3f', Re, n*dt));
        axis equal tight;
        drawnow;
    end
end

fprintf('Simulation finished.\n');

% Optional: save final figure
% saveas(gcf, 'lid_driven_cavity_result.png');
