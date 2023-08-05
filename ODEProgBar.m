function status=ODEProgBar(t,y,flag)
    persistent tf tstart;
    
    if isempty(flag)
        % Integration steps
        ts=mean(t);
        progress=100*ts/tf;
        TextProgressBar(progress);
        status=0;
    else
        switch flag
            case 'init'     % Initializing progress bar
                tstart=tic;
                tf=max(t);
                TextProgressBar('ODE integration: ');
            case 'done'     % Finishing status function
                tf=[];
                TextProgressBar('');
                display([ '   Integration time: ' num2str(toc(tstart))]);
                tstart=[];
            otherwise
                error('odetpbar:UnknownError',...
                    'Unknown error has occured');
        end
    end
end