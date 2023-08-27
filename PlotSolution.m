function PlotSolution(dims,U,x,T,ui)

if (dims == 1)
    hold on
    plot(x,U(1,ui),'linewidth',2);
    plot(x,U(end/10,ui),'linewidth',2);
    plot(x,U(end/2,ui),'linewidth',2)
    plot(x,U(end,ui),'linewidth',2)

    h = legend('$t=0$','$t=10\%$','$t=50\%$','$t=100\%$','interpreter','latex');

    set(gca,'fontsize',24);

elseif (dims == 2)
    M = 6;
    Is = [(1:M-1)*round(length(T)/M), size(U,1)];
    m = round(sqrt(numel(ui)));
    clims = [min(U(:,ui),[],'all'),max(U(:,ui),[],'all')];

    for i = Is
            nexttile();
            imagesc(reshape(U(i,ui),m,m));
            axis equal
            axis tight
            title(['$t=', num2str(T(i)),'$'],'interpreter','latex');
            c = colorbar;
            c.TickLabelInterpreter = 'latex';
            caxis(clims);
    end

end

end