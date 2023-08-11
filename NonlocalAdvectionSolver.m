% This code solves the nonlocal advection system on a periodic 1D domain

if(~exist('setup','var'))
    clear;
    SetupBaseParams;
    % Domain length.
    L = 50;

    % Parameters in the reaction kinetics.
    a = 1; b = 0.4; c=0.5;

    % Nonlocal interaction strength and range
    mu = 20;
    xi = 1;

    % Time interval to solve the equations on.
    T = linspace(0,300,1000);


end
% Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% Spatial step size.
dx = L/(m-1);

% (Sparse) Laplacian and Divergence (Advcection) matrix.
e = ones(m,1);
Lap = spdiags([e,-2*e,e],[1,0,-1],m,m);
Adv = spdiags([e,-e],[1,-1],m,m);

% Periodic boundary conditions.
Lap(1,end) = 1; Lap(end,1) = 1;
Lap = (1/dx)^2*Lap;

Adv(end,1) = 1; Adv(1,end) = -1;
Adv = (1/(2*dx))*Adv;

% Interaction kernel

% the middle element is 0, then with increasing/decreasing increments of
% dx on either side
kernel = cat(2, flip(-dx:-dx:(-L/2 + dx)),0:dx:(L/2 - dx));
% exponentially decaying interaction kernel
kernel = sign(kernel).*exp(-abs(kernel)/xi);
% flip kernel as it get flipped in the convolution
kernel = flip(kernel);
kernel = transpose(kernel);

% Indices corresponding to u variable. This is for plotting code to work.
ui = 1:N;

% Put together the diffusion + kinetics + nonlocal advection terms into a big vector
F = @(t,U) Lap*U + a*U.*(1 - U/c).*(U - b) - dx*(mu/xi)*Adv*(U.*(1 - U).*fft_convolve(U,kernel));

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = c + 1e-2*randn(N,1);

% Setting tolerances for ODE solver
odeOptions = odeset('RelTol',tolerance,'AbsTol',tolerance,'InitialStep',1e-6);
if (showProgressBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[T,U] = ode15s(F,T,U0,odeOptions);

% Function for calculating the circular convolution using FFT
function f = fft_convolve(u, k)
    f = circshift(real(ifft( fft(u).*fft(k, size(u,1)))),-floorDiv(size(k,1)-1,2));
end