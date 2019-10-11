lower = sdpvar(1,1);
rng( 40, 'twister');
N = 64;
b = N/4;
x = sdpvar(N, 1);

f_rosenbrock  = sum(100*(x(2:end ) - x(1:end-1) .^2).^2 +  (1 -x(1:end-1 )).^2);


mBasis = monolist(x(1:b),2);
%Q = rand(length(mBasis));

%Q = hadamard(
%c = rand(length(mBasis) ,1);
%Q = Q*Q';
%f_P = (Q'*mBasis).^2  + c'*mBasis;
%f_P = (Q'*mBasis).^2 ;
Q_P = gallery('lehmer', length(mBasis));
f_P = mBasis' *Q_P * mBasis;

f = f_rosenbrock + f_P;
 
F = sos(f - lower);
obj =  -lower;

opts_csp = sdpsettings('solver','SEDUMI');
%[model, recoverymodel] = export(F , obj, opts);

opts_csp = opts;
opts_csp.sos.csp = 1;
[model_csp, recoverymodel_csp] = export(F , obj, opts_csp);

[x_csp, y_csp, info_csp] = sedumi(model_csp.A', model_csp.b, model_csp.C, model_csp.K);
cost_csp = model_csp.C'*x_csp;
save('polynomial_big_clique.mat', 'F', 'obj', 'x', 'b', 'model_csp')