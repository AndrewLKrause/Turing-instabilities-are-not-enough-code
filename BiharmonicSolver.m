clear;
%This code solves the reaction-diffusion system on a square
%domain in 1D and 2D.

%Seed the random number generator
rng('default');

%Set the spatial dimension to be 1D or 2D.
dimensions=2;

%Number of gridpoints per dimension. Use 60-300 or so for 2D, and ideally
%300-3000 or so for 1D depending on the structures that need to be
%resolved.
m = 100;

%Total number of gridpoints; varies by dimension as:
%2D make N=m^2; 1D make N=m;
if(dimensions==1)
    N = m;
elseif(dimensions==2)
    N = m^2;
end

%Domain length
L =50;

%Parameters in the reaction kinetics
a = 5; b = 0.9;

%Diffusion coefficients
D = 1.55;

%Spatial step size
dx = L/(m-1);

%Time interval to solve the equations on
T = linspace(0,300,1000);

%Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% (Sparse) Laplacian matrix
e = ones(m,1);
Lap = spdiags([e,-2*e,e],[1,0,-1],m,m);


% Neumann boundary conditions
%Lap(1,1) = -1; Lap(end,end) = -1;
% Periodic boundary conditions
Lap(1,end) = 1; Lap(end,1) = 1;

if(dimensions==1)
    %1D Laplacian
    Lap = (1/dx)^2*Lap;
elseif(dimensions==2)
    %2D Laplacian
    I = speye(m);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));
end

%Define the biharmonic operator.
Bih = Lap*Lap;

%Indices corresponding to u variable. This is for plotting code to work.
ui = 1:N; 


%Put together the reaction kinetics+diffusion terms into a big vector
F = @(t,U)-D*Lap*U-Bih*U+a*U.*(1-U).*(U-b);

%Initial condition - this is a small normally distributed perturbation of
%the homogeneous steady state of our kinetics
U0 = 1+1e-3*randn(N,1);

%This is the Jacobian sparsity pattern. That is, if you compute the
%Jacobian of the vector function F above for the vector argument U, this
%matrix is where all of the nonzero elements are. This is important for
%implicit timestepping!
JacSparse = sparse(Bih);
odeOptions = odeset('JPattern',JacSparse,'RelTol',1e-9,'AbsTol',1e-9);

%Solve the system using an implicit stiff timestepper.
[T,U] = ode15s(F,T,U0,odeOptions);
U = U';

