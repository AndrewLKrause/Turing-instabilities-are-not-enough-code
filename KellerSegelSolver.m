function [U,x,ui,vi,uss] = ...
    KellerSegelSolver(dims, m, params, tols, T,showProgBar)
% This code solves the Keller-Segel chemotaxis system on a square
% domain in 1D and 2D.

% Parameters of the model.
[L, a, b, c, d, D] = deal(params{:});

% Spatial step size
dx = L/(m-1);

% Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% (Sparse) Laplacian matrix
eVec = ones(m,1);
Lap = spdiags([eVec,-2*eVec,eVec],[1,0,-1],m,m);
Adv = spdiags([eVec,-eVec],[1,-1],m,m);

% Neumann boundary conditions
% Lap(1,1) = -1; Lap(end,end) = -1;
% Adv(1,2)=0; Adv(end,end-1)=0;
% Periodic boundary conditions
Lap(1,end) = 1; Lap(end,1) = 1;
Adv(end,1) = 1; Adv(1,end) = -1;

if(dims==1)
    N = m;
elseif(dims==2)
    N=m^2;
end
% Indices corresponding to u variable and v variable. THis lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N;

% Reaction kinetics
f = @(u,v) u.*(b - u).*(u - d);
g = @(u,v) u - a*v;

% Set the Laplacian and the nonlinear diffusion DivGrad based on dimension.
% Write the full RHS based on dimension.
if (dims == 1)
    % 1D Laplacian
    Adv = (1/(2*dx))*Adv;
    Lap = (1/dx)^2*Lap;
    
    F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui) -...
        c*(U(ui).*(Lap*U(vi)) + (Adv*U(ui)).*(Adv*U(vi)));
        g(U(ui),U(vi)) + D*Lap*U(vi)];

elseif (dims == 2)
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

%The u value of the HSS used to check for pattern formation
uss = b;

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
JacSparse = sparse([Lap, Lap; speye(N), Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(F,T,U0,odeOptions);
end