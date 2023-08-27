function [U, x, T, ui, vi] = RunAndPlot(modelName,dims)
% Run and plot the model given by the string.

    [U,x,T,ui,vi] = RunSim(modelName,dims);

    PlotSim(dims,U,x,T,ui,vi);

end
