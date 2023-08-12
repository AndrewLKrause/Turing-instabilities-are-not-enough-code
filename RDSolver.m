function [U,x,ui,vi,uss] = ...
    RDSolver(dims, m, params, tols, T,showProgBar)
% This code solves the reaction-diffusion system on a square
% domain in 1D and 2D.

% Parameters of the model.
[L, a, b, c, d, e, D] = deal(params{:});

% Spatial step size
dx = L/(m-1);

% Spatial domain (needed for plotting only).
x = linspace(0,L,m);

% (Sparse) Laplacian matrix.
eVec = ones(m,1);
Lap = spdiags([eVec,-2*eVec,eVec],[1,0,-1],m,m);

% Neumann boundary conditions.
% Lap(1,1) = -1; Lap(end,end) = -1;
% Periodic boundary conditions.
Lap(1,end) = 1; Lap(end,1) = 1;

% Total number of gridpoints; varies by dimension as does Laplacian.
if (dims == 1)
    % 1D Laplacian
    Lap = (1/dx)^2*Lap;
    N = m;
elseif (dims == 2)
    % 2D Laplacian
    I = speye(m);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));
    N = m^2;
end

% Indices corresponding to u variable and v variable. This lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N;

% Reaction kinetics
f = @(u,v) u - v - e*u.^3;
g = @(u,v) a*v.*(v + c).*(v - d) + b*u - e*v.^3;

% Put together the reaction kinetics+diffusion terms into a big vector.
F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui);
    g(U(ui),U(vi)) + D*Lap*U(vi)];

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics.
U0 = 1e-2*randn(2*N,1);

%The u value of the HSS used to check for pattern formation
uss = 0;

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
JacSparse = sparse([Lap,speye(N);speye(N),Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(F,T,U0,odeOptions);
end