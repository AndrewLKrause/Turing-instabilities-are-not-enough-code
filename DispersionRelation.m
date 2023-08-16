function lambda = DispersionRelation(modelName, params)
% This function determines if the Turing Conditions are satisfied for the
% given model and parameters.

switch modelName
    case 'RD'
        [~, a, b, c, d, ~, D] = deal(params{:});
        fu = 1; fv = -1; gu = b; gv = -a*c*d;

        trM = fu+gv-rho_k*(1+D);
        detM = rho_k^2*D-rho_k*(gv+D*fu)+fu*gv-fv*gu;
        lambda = (trM+sqrt(trM^2-4*detM))/2;
        %rho_k^2*(d1*d4-d2*d3)+rho_k*(d2*gu+d3*fv-d1*gv-d4*fu))+fu*gv-fv*gu;
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
        [L, a, b, c,d, D] = deal(params{:});
        lambda = @(k)-a*c*(c-b)-D*k^2 ...
            + 2*pi*c*(1-c)*d*(1/(2*pi))*(k^2)/((1 + k^2)^(3/2));
        TuringUnstable = max(arrayfun(lambda,(0:1000)*pi/L))>0;
    otherwise
        disp('Unknown method.')
        return;
end


end