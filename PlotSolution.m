function PlotSolution(dims,U, x, ui)

close all;
if(dims==1)
    plot(x,U(1,ui),'linewidth',2); hold on
    plot(x,U(end/10,ui),'linewidth',2);
    plot(x,U(end/2,ui),'linewidth',2)
    plot(x,U(end,ui),'linewidth',2)

    h = legend('$t=0$','$t=10\%$','$t=50\%$','$t=100\%$','interpreter','latex');

    set(gca,'fontsize',24);

elseif(dims==2)
    M=6;
    Is = round(length(T)/M);
    for i = 1:M-1
            figure; imagesc(reshape(U(Is*i,ui),m,m)); title(['$t=', num2str(T(Is*i)),'$'],'interpreter','latex'); colorbar;
    end
    figure; imagesc(reshape(U(end,ui),m,m)); title(['$t=', num2str(T(end)),'$'],'interpreter','latex'); colorbar;

end

end