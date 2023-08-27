function PlotKymograph(U,x,T,ui)
% Plots a space-time diagram (Kymograph) of u.

imagesc(T,x,(U(:,ui))'); 
set(gca,'YDir','normal')
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
colormap(viridis)
shading interp

c = colorbar;
c.Label.String = '$u$';
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';

set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24);

end