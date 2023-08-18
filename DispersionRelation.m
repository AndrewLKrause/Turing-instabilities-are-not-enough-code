function lambda = DispersionRelation(modelName, params,rho_k)
% This function determines if the Turing Conditions are satisfied for the
% given model and parameters.

switch modelName
    case 'RD'
        [~, a, b, c, d, ~, D] = deal(params{:});
        fu = 1; fv = -1; gu = b; gv = -a*c*d;
        trM = fu+gv-rho_k*(1+D);
        detM = rho_k^2*D-rho_k*(gv+D*fu)+fu*gv-fv*gu;
        lambda = (trM+sqrt(trM^2-4*detM))/2;
    case 'KellerSegel'
        [~, a, b, c, d, D] = deal(params{:});
        fu = -b*(b-d); fv = 0; gu = 1; gv = -a;
        trM = fu+gv-rho_k*(1+D);
        detM = rho_k^2*D+-rho_k*(c*b*gu+gv+D*fu)+fu*gv-fv*gu;
        lambda = (trM+sqrt(trM^2-4*detM))/2;
    case 'Biharmonic'
        [~, a, b, c, D] = deal(params{:});
        lambda = D*rho_k-rho_k^2+a*c*(b-c);
    case 'NonlocalAdvection'
        [~, a, b, c,d, D] = deal(params{:});
        lambda = -a*c*(c-b)-D*rho_k + c*(1-c)*d*(rho_k)/((1 + rho_k)^(3/2));
    otherwise
        disp('Unknown method.')
        return;
end


end