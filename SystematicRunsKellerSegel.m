% Set the spatial dimension to be 1D or 2D.
dimensions = 1;

% Show a progress bar?
showProgressBar = false;

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

setup = 1;

NumRuns = 1000; % Number of runs to check.
Var = 0.05; R = @()(1-Var+2*Var*rand);

% Time interval to solve the equations on.
T = linspace(0,10000,10);
Patterning = zeros(NumRuns,1);
Patterns = zeros(NumRuns,1);
if (~showProgressBar)
    TextProgressBar('Running: ');
end
sims = cell(NumRuns,1);

for iRun = 1:NumRuns
    % Seed the random number generator.
    rng(iRun);

    % Domain length
    L = 80*R();

    % Parameters in the reaction kinetics
    c = 3*R(); A = 0.8*R(); a = 1*R();

    % Diffusion coefficients
    Du = 1*R();
    Dv = 1*R();

    % Spatial step size.
    dx = L/(m-1);

    KellerSegelSolver;
    if (~showProgressBar)
        TextProgressBar(100 * iRun / NumRuns)
    end
    Patterns(iRun) = max(U(end,ui))-min(U(end,ui));
    Patterning(iRun) = max(U(end,ui));

    % Save.
    runDetails = struct();
    runDetails.L = L;
    runDetails.c = c;
    runDetails.a = a;
    runDetails.A = A;
    runDetails.Du = Du;
    runDetails.Dv = Dv;
    runDetails.dx = dx;
    runDetails.U0 = U0;
    runDetails.Patterns = Patterns(iRun);
    runDetails.Patterning = Patterning(iRun);
    sims{iRun} = runDetails;
end
TextProgressBar('')