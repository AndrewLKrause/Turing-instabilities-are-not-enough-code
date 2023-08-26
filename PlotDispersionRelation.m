function PlotDispersionRelation(modelName,max_k)
% Plot a dispersion relation corresponding to a model with the standard
% parameters given in CreateBaseParams.

[~,~, params, ~] = CreateBaseParams(modelName, 1);

L = params{1};

rho = linspace(0,(max_k*pi/L)^2,1e4);
lambda = DispersionRelation(modelName, params,rho);

close all;

plot(rho,lambda,'linewidth',2);hold on;
xlabel('$\rho_k$','interpreter','latex')
ylabel('$\lambda_k$','interpreter','latex')

rho_k = ((0:max_k)*(pi/L)).^2;
lambda_k = DispersionRelation(modelName, params,rho_k)

line([0, (max_k*pi/L)^2], [0, 0], 'color','k');
plot(rho_k,lambda_k,'.','markersize',15);
axis tight;
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24);

end