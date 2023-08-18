function [U, x, T, ui, vi] = RunAndPlot(modelName,dims)
% Run and plot the model given by the string.

% Show a progress bar?
showProgBar = true;
% Set default random seed and get numerical parameters.
rng('default');
[m,tols, BaseParams, Solver] = CreateBaseParams(modelName, dims);

switch modelName
    case 'RD'
        T = linspace(0,330,1000);% Solution timescale.
    case 'KellerSegel'
        T = linspace(0,50,1000);
    case 'Biharmonic'
        T = linspace(0,180,1000);
    case 'NonlocalAdvection'
        T = linspace(0,18,1000);
end

[U,x,ui,vi] = Solver(dims, m, BaseParams, tols, T,showProgBar);

if(dims==1)
    PlotKymograph(U,x,T,ui);
elseif(dims==2)
    PlotSolution(dims,U, x, ui)
end
%AnimateSolution(dims,U, x,T, ui,vi)
