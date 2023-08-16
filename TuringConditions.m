function TuringUnstable = TuringConditions(modelName, params)
% This function determines if the Turing Conditions are satisfied for the
% given model and parameters. Note this incorporates domain length/finite 
% wavemode selection by looking at the first 500 spatial eigenvalues.

max_k = 500;
rho_k = ((0:max_k)*(pi/L)).^2;

lambda = DispersionRelation(modelName, params,rho_k);

TuringUnstable = (real(lambda(1)) < 0) && (max(real(lambda))>0);

end