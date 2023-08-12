function TuringUnstable = TuringConditions(modelName, params)
% This function determines if the Turing Conditions are satisfied for the
% given model and parameters.

switch modelName
    case 'RD'
        [L, a, b, c, d, e, D] = deal(params{:});
        fu = 1; fv = -1; gu = b; gv = -a*c*d;
        TuringUnstable = (fu+gv<0) && (fu*gv-fv*gu>0) && (D*fu+gv>0) ...
            && ((D*fu+gv)^2-4*D*(fu*gv-gu*fv)>0);
    case 'KellerSegel'
        [L, a, b, c, d, D] = deal(params{:});
        fu = b*(b-d); fv = 0; gu = 1; gv = -a;
        TuringUnstable = (D*fu+gv+c*b*gu>0)...
            && ((D*fu+gv+c*b*gu)^2-4*D*(fu*gv-gu*fv)>0);
    case 'Biharmonic'
        [L, a, b, c, D] = deal(params{:});
        TuringUnstable = D^2-4*a*c*(b-c)>0;
    case 'NonlocalAdvection'

    otherwise
        disp('Unknown method.')
        return;
end


end