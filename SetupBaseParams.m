% Seed the random number generator.
rng('default');

% Set the spatial dimension to be 1D or 2D.
dimensions = 1;

% Show a progress bar?
showProgressBar = true;

% Numerical tolerances (absolute and relative).
tolerance = 1e-9;

% Number of gridpoints per dimension. Use 60-300 or so for 2D, and ideally
% 300-3000 or so for 1D depending on the structures that need to be
% resolved.
m = 1000;

% Total number of gridpoints; varies by dimension as:
% 2D make N=m^2; 1D make N=m;
if (dimensions == 1)
    N = m;
elseif (dimensions == 2)
    N = m^2;
end