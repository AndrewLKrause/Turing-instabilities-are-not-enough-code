clear;
% This code solves the reaction-diffusion system on a square
% domain in 1D and 2D.

% Seed the random number generator
rng('default');

% Set the spatial dimension to be 1D or 2D.
dimensions = 2;

% Show a progress bar?
showProgressBar = true;

% Numerical tolerances (absolute and relative).
tolerance = 1e-9;

% Number of gridpoints per dimension. Use 60-300 or so for 2D, and ideally
% 300-3000 or so for 1D depending on the structures that need to be
% resolved.
m = 200;

% Total number of gridpoints; varies by dimension as:
% 2D make N=m^2; 1D make N=m.
if (dimensions == 1)
    N = m;
elseif (dimensions == 2)
    N = m^2;
end

% Domain length
L = 100;

% Parameters in the reaction kinetics
epsilon = 0.01;
a = 7/4; b = 10; c = (7972/4067);

% Diffusion coefficients
Du = 1;
Dv = 25;

% Spatial step size
dx = L/(m-1);

% Time interval to solve the equations on
T = linspace(0,1000,1000);

% Spatial domain (needed for plotting only).
x = linspace(0,L,m);

% (Sparse) Laplacian matrix.
e = ones(m,1);
Lap = spdiags([e,-2*e,e],[1,0,-1],m,m);

% Neumann boundary conditions.
% Lap(1,1) = -1; Lap(end,end) = -1;
% Periodic boundary conditions.
Lap(1,end) = 1; Lap(end,1) = 1;

if (dimensions == 1)
    % 1D Laplacian
    Lap = (1/dx)^2*Lap;
elseif (dimensions == 2)
    % 2D Laplacian
    I = speye(m);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));
end

% Indices corresponding to u variable and v variable. This lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N;

% Reaction kinetics
f = @(u,v) u - v - epsilon*u.^3;
g = @(u,v) a*v.*(v + c).*(v - 5) + a*b*u - epsilon*v.^3;

% Put together the reaction kinetics+diffusion terms into a big vector.
F = @(t,U)[f(U(ui),U(vi)) + Du*Lap*U(ui);
           g(U(ui),U(vi)) + Dv*Lap*U(vi)];

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics.
U0 = 1e-3*randn(2*N,1);

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
JacSparse = sparse([Lap,speye(N);speye(N),Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tolerance,'AbsTol',tolerance);
if (showProgressBar) 
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[T,U] = ode15s(F,T,U0,odeOptions);
