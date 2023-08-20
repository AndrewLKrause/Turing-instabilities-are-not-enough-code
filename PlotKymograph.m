function PlotKymograph(U,x,T,ui)
% Plots a space-time diagram (Kymograph) of u.


close all;

imagesc(T,x,(U(:,ui))'); set(gca,'YDir','normal')%ax.YTickLabel = flip(ax.YTickLabel);
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
colormap(viridis)

c = colorbar; c.TickLabelInterpreter='latex';

set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24);

end