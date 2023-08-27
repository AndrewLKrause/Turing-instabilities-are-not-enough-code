function [TuringUnstable, k] = TuringConditions(modelName, params, dims)
% This function determines if the Turing Conditions are satisfied for the
% given model and parameters. Note this incorporates domain length/finite 
% wavemode selection by looking at the first 500 spatial eigenvalues.

if nargin < 3
    dims = 1;
end

L = params{1};
max_k = 500;
ks = 0:max_k;
rho_k = (ks*(pi/L)).^2;

lambda = DispersionRelation(modelName,params,rho_k,dims);

TuringUnstable = (real(lambda(1)) < 0) && (max(real(lambda))>0);
k = NaN;
if TuringUnstable
    k = ks(find(real(lambda) > 0,1,'first'));
end

end