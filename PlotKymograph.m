function PlotKymograph(U,x,T,ui)
% Plots a space-time diagram (Kymograph) of u.


close all;

imagesc(x,T,(U(:,ui))); set(gca,'YDir','normal')%ax.YTickLabel = flip(ax.YTickLabel);
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
colorbar;

set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24);

end