function PlotDispersionRelation(modelName,max_k,dims)
% Plot a dispersion relation corresponding to a model with the standard
% parameters given in CreateBaseParams.
    if nargin < 3
        dims = 1;
    end

    [~,~, params, ~] = CreateBaseParams(modelName, dims);

    L = params{1};

    rho = linspace(0,(max_k*pi/L)^2,1e4);
    lambda = DispersionRelation(modelName, params, rho, dims);

    plot(rho,lambda,'linewidth',2);
    hold on;
    xlabel('$\rho_k$','interpreter','latex')
    ylabel('$\lambda_k$','interpreter','latex')

    rho_k = ((0:max_k)*(pi/L)).^2;
    lambda_k = DispersionRelation(modelName, params, rho_k, dims);

    line([0, (max_k*pi/L)^2], [0, 0], 'color','k');
    plot(rho_k,lambda_k,'.','markersize',10, 'color','#D95319');
    axis tight;
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'fontsize',24);

end