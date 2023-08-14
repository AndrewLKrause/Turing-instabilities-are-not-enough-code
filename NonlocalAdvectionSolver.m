function [U,x,ui,vi,uss] = ...
    NonlocalAdvectionSolver(dims, m, params, tols, T,showProgBar)
% This code solves the nonlocal advection system on a periodic 1D domain.

% Parameters of the model.
[L, a, b, c, d, D] = deal(params{:});

% Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% Spatial step size.
dx = L/(m-1);

% (Sparse) Laplacian and Divergence (Advcection) matrix.
eVec = ones(m,1);
Lap = spdiags([eVec,-2*eVec,eVec],[1,0,-1],m,m);
Adv = spdiags([eVec,-eVec],[1,-1],m,m);

% Periodic boundary conditions.
Lap(1,end) = 1; Lap(end,1) = 1;
Lap = (1/dx)^2*Lap;

Adv(end,1) = 1; Adv(1,end) = -1;
Adv = (1/(2*dx))*Adv;

% Interaction kernel the middle element is 0, then with 
% increasing/decreasing increments of dx on either side.
kernel = cat(2, flip(-dx:-dx:(-L/2 + dx)),0:dx:(L/2 - dx));
% Exponentially decaying interaction kernel.
kernel = sign(kernel).*exp(-abs(kernel));
% Flip kernel as it get flipped in the convolution.
kernel = flip(kernel);
kernel = transpose(kernel);

% Function for calculating the circular convolution using FFT
fft_convolve = @(u)circshift(real(ifft(fft(u).*fft(kernel, size(u,1)))),-floor(size(kernel,1)-1)/2);

% Indices corresponding to u variable. This is for plotting code to work.
ui = 1:m; vi=[];

% Put together the diffusion + kinetics + nonlocal advection terms.
F = @(t,U) D*Lap*U + a*U.*(1 - U/c).*(U - b) ...
    - dx*(d)*Adv*(U.*(1 - U).*fft_convolve(U));

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics.
U0 = c + 1e-2*randn(m,1);

%The u value of the HSS used to check for pattern formation.
uss = b;

% Setting tolerances for ODE solver.
odeOptions = odeset('RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(F,T,U0,odeOptions);

end