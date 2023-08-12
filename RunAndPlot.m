function RunAndPlot(modelName,dims)
% Run and plot the model given by the string.

% Show a progress bar?
showProgBar = true;
% Set default random seed and get numerical parameters.
rng('default');
[m,tols] = CreateBaseParams(dims);

switch modelName
    case 'RD'
        BaseParams = {100, 1.75, 18, 2, 5, 0.02,25};
        %        {  L,   a,  b, c, d,    e, D}
        Solver = @RDSolver;
        T = linspace(0,300,1000);% Solution timescale.
    case 'KellerSegel'
        BaseParams = {80, 1, 1, 3, 0.8,1};
        %            { L, a, b, c,   d, D}
        Solver = @KellerSegelSolver;
        T = linspace(0,1000,1000);
    case 'Biharmonic'
        BaseParams = {100, 5, 0.9, 1, 1.45};
        %            {  L, a,   b, c,    D}
        Solver = @BiharmonicSolver;
        T = linspace(0,200,1000);
    case 'NonlocalAdvection'
        if(dims~=1)
            disp('Currently only 1D is implemented for this model.');
            return;
        end
        BaseParams = {30, 1, 0.45, 0.5, 20, 1};
        %            { L, a,    b,   c,  d, D}
        Solver = @NonlocalAdvectionSolver;
        T = linspace(0,25,1000);
    otherwise
        disp('Unknown method.')
        return;
end

[U,x,ui,vi] = Solver(dims, m, BaseParams, tols, T,showProgBar);
%PlotSolution(dims,U, x, ui)
AnimateSolution(dims,U, x,T, ui,vi)
