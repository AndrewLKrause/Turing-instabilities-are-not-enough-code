clear;
% This code solves the reaction-diffusion system on a square
% domain in 1D and 2D.

% Seed the random number generator.
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
% 2D make N=m^2; 1D make N=m;
if (dimensions == 1)
    N = m;
elseif (dimensions == 2)
    N = m^2;
end

% Domain length
L = 80;

% Parameters in the reaction kinetics
c = 3; A = 0.5; a = 1;

% Diffusion coefficients
Du = 1;
Dv = 1;

% Spatial step size
dx = L/(m-1);

% Time interval to solve the equations on
T = linspace(0,500,1000);

% Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% (Sparse) Laplacian matrix
e = ones(m,1);
Lap = spdiags([e,-2*e,e],[1,0,-1],m,m);
Adv = spdiags([e,-e],[1,-1],m,m);

% Neumann boundary conditions
% Lap(1,1) = -1; Lap(end,end) = -1;
% Adv(1,2)=0; Adv(end,end-1)=0;
% Periodic boundary conditions
Lap(1,end) = 1; Lap(end,1) = 1;
Adv(end,1) = 1; Adv(1,end) = -1;


% Indices corresponding to u variable and v variable. THis lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N;

% Reaction kinetics
f = @(u,v) u.*(1 - u).*(u - A);
g = @(u,v) u - a*v;

% Set the Laplacian and the nonlinear diffusion DivGrad based on dimension.
% Write the full RHS based on dimension.
if (dimensions == 1)
    % 1D Laplacian
    Adv = (1/(2*dx))*Adv;
    Lap = (1/dx)^2*Lap;
    
    F = @(t,U)[f(U(ui),U(vi)) + Du*Lap*U(ui) -...
               c*(U(ui).*(Lap*U(vi)) + (Adv*U(ui)).*(Adv*U(vi)));
               g(U(ui),U(vi)) + Dv*Lap*U(vi)];
    
elseif (dimensions == 2)
    % 2D Laplacian
    I = speye(m);
    Advx = (1/(2*dx))*kron(Adv,I);
    Advy = (1/(2*dx))*kron(I, Adv);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));

    F = @(t,U)[f(U(ui),U(vi)) + Du*Lap*U(ui) -...
               c*(U(ui).*(Lap*U(vi)) +...
                (Advx*U(ui)).*(Advx*U(vi)) + (Advy*U(ui)).*(Advy*U(vi)));
               g(U(ui),U(vi)) + Dv*Lap*U(vi)];

end


% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = [1 + 1e-3*randn(N,1); 1/a + 1e-3*randn(N,1)];

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
JacSparse = sparse([Lap, Lap; speye(N), Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tolerance,'AbsTol',tolerance);
if (showProgressBar) 
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[T,U] = ode15s(F,T,U0,odeOptions);
