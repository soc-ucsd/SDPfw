load('LR_120.mat')

model = model_cons_trans;

[F, obj] = sedumi2yalmip(model.A',model.b,model.c,model.K);

opts_csp = sdpsettings('solver','sedumi', 'verbose', 0);
opts_csp.sos.csp = 1;
opts_csp.sos.model = 2;

%Dualize it: software for automatic primal and dual conversions of
%conic programsDualize it: software for automatic primal and dual conversions of
%conic programs

%dualize (primalize?)
[Fd,objd,X,t,err] = dualize(F,obj,0);
[model_dual, ~] = export(Fd, objd, opts_csp);

%model_dual_unc = model_dual;
%[x_dual, y_dual, info_dual] = sedumi(-model_dual.A, model_dual.b,model_dual.C,  model_dual.K);
%res_dual_unc = convert_sedumi2mosek(-model_dual.A',model_dual.b,model_dual.C,model_dual.K);

%save('LR_120.mat', '-append', 'res_dual_unc', 'model_dual_unc')