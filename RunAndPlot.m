function RunAndPlot(modelName,dims)
% Run and plot the model given by the string.

% Show a progress bar?
showProgBar = true;
% Set default random seed and get numerical parameters.
rng('default');
[m,tols, BaseParams, Solver] = CreateBaseParams(modelName, dims);

switch modelName
    case 'RD'
        T = linspace(0,300,1000);% Solution timescale.
    case 'KellerSegel'
        T = linspace(0,1000,1000);
    case 'Biharmonic'
        T = linspace(0,200,1000);
    case 'NonlocalAdvection'
        T = linspace(0,25,1000);
end

[U,x,ui,vi] = Solver(dims, m, BaseParams, tols, T,showProgBar);
%PlotSolution(dims,U, x, ui)
AnimateSolution(dims,U, x,T, ui,vi)
