function [U,x,T,ui,vi] = RunSim(modelName,dims)

    % Show a progress bar?
    showProgBar = true;
    % Set default random seed and get numerical parameters.
    rng('default');
    [m,tols, BaseParams, Solver] = CreateBaseParams(modelName, dims);

    switch modelName
        case 'RD'
            BaseParams = {100, 1.75, 18, 2, 5, 0.02,30};
            %        {  L,   a,  b, c, d,    e, D}
            Solver = @RDSolver;
            T = linspace(0,330,1000);% Solution timescale.
        case 'KellerSegel'
            BaseParams = {80, 1, 1, 3, 0.8,1};
            %            { L, a, b, c,   d, D}
            Solver = @KellerSegelSolver;
            T = linspace(0,50,1000);
        case 'Biharmonic'
            BaseParams = {100, 5, 0.9, 1, 1.45};
            %            {  L, a,   b, c,    D}
            Solver = @BiharmonicSolver;
            T = linspace(0,200,1000);
        case 'NonlocalAdvection'
            T = linspace(0,18,1000);
    end

    [U,x,ui,vi] = Solver(dims, m, BaseParams, tols, T,showProgBar);

end