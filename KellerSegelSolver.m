% This code solves the Keller-Segel chemotaxis system on a square
% domain in 1D and 2D.

if(~exist('setup','var'))
    clear;
    SetupBaseParams;

    % Domain length
    L = 80;

    % Parameters in the reaction kinetics
    a = 1; b = 1; c = 3; d = 0.8; 

    % Diffusion coefficients
    D = 1;

    % Time interval to solve the equations on
    T = linspace(0,500,1000);
end


% Spatial step size
dx = L/(m-1);

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
f = @(u,v) u.*(b - u).*(u - d);
g = @(u,v) u - a*v;

% Set the Laplacian and the nonlinear diffusion DivGrad based on dimension.
% Write the full RHS based on dimension.
if (dimensions == 1)
    % 1D Laplacian
    Adv = (1/(2*dx))*Adv;
    Lap = (1/dx)^2*Lap;

    F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui) -...
        c*(U(ui).*(Lap*U(vi)) + (Adv*U(ui)).*(Adv*U(vi)));
        g(U(ui),U(vi)) + D*Lap*U(vi)];

elseif (dimensions == 2)
    % 2D Laplacian
    I = speye(m);
    Advx = (1/(2*dx))*kron(Adv,I);
    Advy = (1/(2*dx))*kron(I, Adv);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));

    F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui) -...
        c*(U(ui).*(Lap*U(vi)) +...
        (Advx*U(ui)).*(Advx*U(vi)) + (Advy*U(ui)).*(Advy*U(vi)));
        g(U(ui),U(vi)) + D*Lap*U(vi)];

end


% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = [b + 1e-2*randn(N,1); 1/a + 1e-2*randn(N,1)];

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
JacSparse = sparse([Lap, Lap; speye(N), Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tolerance,'AbsTol',tolerance,'InitialStep',1e-6);
if (showProgressBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[T,U] = ode15s(F,T,U0,odeOptions);
