function [NumPatterns, NumPatterning, NumTuring] = SystematicRuns(modelName,dims)
% This function performs NumRuns simulations of a given model and dimension
% The outputs are the number of sims that had a 'pattern' at the final 
% time, and the number which left the HSS, respectively.

if nargin < 1
    modelName = 'RD';
end

if nargin < 2
    dims = 1;
end

modelName = convertStringsToChars(modelName);

% Show a progress bar?
showProgBar = false;
% Set default random seed and get numerical parameters.
rng('default');
[m,tols, BaseParams, Solver] = CreateBaseParams(modelName, dims);

NumRuns = 100; % Number of runs to check.
Var = 0.05; % Percentage variation from base parameter values

NumParams = length(BaseParams);

% Time interval to solve the model on.
T = linspace(0,10000,10);

Patterning = zeros(NumRuns,1);
Patterns = zeros(NumRuns,1);
TuringUnstable = zeros(NumRuns,1);

if (showProgBar)
    TextProgressBar('Running: ');
end
RunSolutions = cell(NumRuns,1);
RunParameters = cell(NumRuns,1);


% Seed the random number generator and generate a latin hypercube sample.
rng('default');
LHS=(1-Var+2*Var*lhsdesign(NumRuns,NumParams));

%Create a parallel pool or do nothing if one exists.
gcp;

parfor iRun = 1:NumRuns
    % Use a random seed based on the current run.
    rng(iRun);
    params = cell(1,NumParams);
    for jParam = 1:NumParams
        params{jParam} = BaseParams{jParam}*LHS(iRun,jParam);
    end
    % These are parameters (dims,L, a, b, c,D, LHS, d, e)
    [U,~,ui,~,uss] = Solver(dims, m, params, tols, T,false);

    if (showProgBar)
        TextProgressBar(100 * iRun / NumRuns)
    end
    Patterns(iRun) = max(U(end,ui))-min(U(end,ui));
    Patterning(iRun) = max(abs(U(end,ui)-uss))>1e-5; 
    TuringUnstable(iRun) = TuringConditions(modelName, params);
    % Variables to reconstruct runs.
    RunSolutions{iRun} = U;
    RunParameters{iRun} = params;
end

% Output the number of sims that had a 'pattern' at the final time, and the
% number which left the HSS
NumPatterns = nnz(Patterns > 1e-5);
NumPatterning = nnz(Patterning); 
NumTuring = nnz(TuringUnstable);

%NB LHS stores all random numbers except U0.
save(['DataRuns',modelName, num2str(dims), 'D.mat']);


TextProgressBar('')
end