clear;
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

NumRuns = 10000; % Number of runs to check.
Var = 0.05; 

% Time interval to solve the equations on.
T = linspace(0,10000,10);
Patterning = zeros(NumRuns,1);
Patterns = zeros(NumRuns,1);
if (~showProgressBar)
    TextProgressBar('Running: ');
end
sims = cell(NumRuns,1);

% Seed the random number generator and generate a latin hypercube sample.
rng('default');
LHS=(1-Var+2*Var*lhsdesign(NumRuns,8));

for iRun = 1:NumRuns


    % Domain length
    L = 100*LHS(iRun, 1);

    % Parameters in the reaction kinetics
    epsilon = 0.01*LHS(iRun, 2);
    a = 1.75*LHS(iRun, 3); b = 17.5*LHS(iRun, 4); c = 2*LHS(iRun, 5); 
    d = 5*LHS(iRun, 6);

    % Diffusion coefficients
    Du = 1*LHS(iRun, 7);
    Dv = 25*LHS(iRun, 8);

    % Spatial step size.
    dx = L/(m-1);

    RDSolver;
    if (~showProgressBar)
        TextProgressBar(100 * iRun / NumRuns)
    end
    Patterns(iRun) = max(U(end,ui))-min(U(end,ui));
    Patterning(iRun) = max(U(end,ui));

    % Save.
    runDetails = struct();
    runDetails.L = L;
    runDetails.epsilon = epsilon;
    runDetails.a = a;
    runDetails.b = b;
    runDetails.c = c;
    runDetails.Du = Du;
    runDetails.Dv = Dv;
    runDetails.dx = dx;
    runDetails.U0 = U0;
    runDetails.Patterns = Patterns(iRun);
    runDetails.Patterning = Patterning(iRun);
    sims{iRun} = runDetails;
end
save("SystematicRunsRD.mat", "sims");
TextProgressBar('')
