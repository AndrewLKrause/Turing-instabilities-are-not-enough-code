% This code solves the biharmonic system on a square
% domain in 1D and 2D.

if(~exist('setup','var'))
    clear;
    SetupBaseParams;
    % Domain length.
    L = 50;

    % Parameters in the reaction kinetics.
    a = 5; b = 0.9;

    % Diffusion coefficients.
    D = 1.55;

    % Spatial step size.
    dx = L/(m-1);

    % Time interval to solve the equations on.
    T = linspace(0,300,1000);

    % Spatial domain (needed for plotting only)
    x = linspace(0,L,m);
end

% (Sparse) Laplacian matrix.
e = ones(m,1);
Lap = spdiags([e,-2*e,e],[1,0,-1],m,m);

% Neumann boundary conditions.
% Lap(1,1) = -1; Lap(end,end) = -1;
% Periodic boundary conditions.
Lap(1,end) = 1; Lap(end,1) = 1;

if (dimensions == 1)
    % 1D Laplacian.
    Lap = (1/dx)^2*Lap;
elseif (dimensions == 2)
    % 2D Laplacian.
    I = speye(m);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));
end

% Define the biharmonic operator.
Bih = Lap*Lap;

% Indices corresponding to u variable. This is for plotting code to work.
ui = 1:N;

% Put together the reaction kinetics+diffusion terms into a big vector
F = @(t,U) -D*Lap*U - Bih*U + a*U.*(1 - U).*(U - b);

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = 1 + 1e-3*randn(N,1);

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
JacSparse = sparse(Bih);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tolerance,'AbsTol',tolerance);
if (showProgressBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[T,U] = ode15s(F,T,U0,odeOptions);

