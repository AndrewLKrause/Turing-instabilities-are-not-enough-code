close all;
f1=figure;
plotV = exist('vi','var');
if (dimensions == 1)
    maxU = max(max(U)); minU = min(min(U));
    for i=1:length(T)
        if ~ishghandle(f1)
            break
        end
        plot(x,U(i,ui),'linewidth',2); hold on
        if(plotV)
            plot(x,U(i,vi),'--','linewidth',2);
            %legend('$u$','$v$','interpreter','latex');
        end
        set(gca,'fontsize',24);
        axis tight;
        set(gca,'YLim',[minU,maxU]);
        title(['$t = ',num2str(T(i)),'$'])
        pause(5/length(T));
        hold off;
    end
elseif (dimensions==2)
    maxU = max(max(U(:,ui))); minU = min(min(U(:,ui)));
    imagesc(reshape(U(1,ui),m,m)); colorbar; 
    caxis([minU,maxU]); ax = gca; ax.XTick = []; ax.YTick = [];
    hold on;
    for i=1:length(T)
        if ~ishghandle(f1)
            break
        end
        imagesc(reshape(U(i,ui),m,m)); 
        title(['$t = ',num2str(T(i)),'$'])
        pause(5/length(T));
    end

end