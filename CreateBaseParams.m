function [m, tols, BaseParams, Solver] = CreateBaseParams(modelName, dims)

switch modelName
    case 'RD'
        BaseParams = {100, 1.75, 18, 2, 5, 0.02,30};
        %            {  L,   a,  b, c, d,    e, D}
        Solver = @RDSolver;
    case 'KellerSegel'
        BaseParams = {80, 1, 1, 3, 0.8,1};
        %            { L, a, b, c,   d, D}
        Solver = @KellerSegelSolver;
    case 'Biharmonic'
        BaseParams = {100, 5, 0.9, 1, 1.45};
        %        {  L, a,   b, c,    D}
        Solver = @BiharmonicSolver;
    case 'NonlocalAdvection'
        if(dims~=1)
            disp('Only 1D is implemented in MATLAB for this model.');
            return;
        end
        BaseParams = {30, 1, 0.45, 0.5, 20, 1};
        %            { L, a,    b,   c,  d, D}
        Solver = @NonlocalAdvectionSolver;
    otherwise
        disp('Unknown model.')
        return;
end

% Number of gridpoints per dimension. Use 60-300 or so for 2D, and ideally
% 300-3000 or so for 1D depending on the structures that need to be
% resolved.
if(dims==1)
    m = 1000;
elseif(dims==2)
    m=100;
end

% Numerical tolerances (absolute and relative).
tols = 1e-9;

end