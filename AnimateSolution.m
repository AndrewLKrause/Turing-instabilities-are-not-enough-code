function AnimateSolution(dims,U, x,T, ui,vi)
close all;
f1=figure;
plotV = exist('vi','var')&& ~isempty(vi);
if (dims == 1)
    maxU = max(max(U)); minU = min(min(U));
    for i=1:length(T)
        if ~ishghandle(f1)
            break
        end
        hold off
        plot(x,U(i,ui),'linewidth',2); hold on
        if(plotV)
            plot(x,U(i,vi),'--','linewidth',2);
            %legend('$u$','$v$','interpreter','latex');
        end
        set(gca,'fontsize',24);
        axis tight;
        set(gca,'YLim',[minU,maxU]);
        title(['$t = ',num2str(T(i)),'$'],'interpreter','latex')
        pause(2/length(T));
    end
elseif (dims==2)
    m = sqrt(length(ui));
    maxU = max(max(U(:,ui))); minU = min(min(U(:,ui)));
    imagesc(reshape(U(1,ui),m,m)); colorbar; 
    caxis([minU,maxU]); ax = gca; ax.XTick = []; ax.YTick = [];
    for i=1:length(T)
        if ~ishghandle(f1)
            break
        end
        imagesc(reshape(U(i,ui),m,m)); 
        title(['$t = ',num2str(T(i)),'$'],'interpreter','latex')
        pause(1/length(T));
    end
end

end