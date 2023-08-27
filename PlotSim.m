function PlotSim(dims,U,x,T,ui,vi)

    if (dims == 1)
        PlotKymograph(U,x,T,ui);
    elseif (dims == 2)
        PlotSolution(dims,U,x,T,ui)
    end

end