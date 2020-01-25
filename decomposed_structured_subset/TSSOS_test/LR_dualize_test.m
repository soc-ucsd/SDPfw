%load('LR_24.mat')
load('LR_120.mat', 'model_cons_trans')
model = model_cons_trans;
%model_dual = model_cons_trans_dual;

[F, obj] = sedumi2yalmip(model.A',model.b,model.c,model.K);

opts = sdpsettings('solver','sedumi', 'verbose', 0);
opts.sos.csp = 1;
opts.sos.model = 2;

%Dualize it: software for automatic primal and dual conversions of
%conic programsDualize it: software for automatic primal and dual conversions of
%conic programs

%dualize (primalize?)
[Fd,objd,X,t,err] = dualize(F,obj,0);
[model_dual_cons_trans, ~] = export(Fd, objd, opts);

%maybe?
model_dual_cons_trans.A = -model_dual_cons_trans.A;


%model_dual_unc = model_dual;
%CSP
%[x, y, info] = sedumi(model.A', model.b, model.c, model.K);
%cost = -model.c'*x;

% prob = convert_sedumi2mosek(model.A', model.b, model.c, model.K);
% [r, res] = mosekopt('minimize', prob);
% 
% 
% prob_dual = convert_sedumi2mosek(-model_dual.A, model_dual.b,model_dual.C,  model_dual.K);
% [r_dual, res_dual] = mosekopt('minimize', prob_dual);
% 

%Dualized
%[x_dual, y_dual, info_dual] = sedumi(-model_dual.A, model_dual.b,model_dual.C,  model_dual.K);
%cost_dual = model_dual.C'*x_dual;

%[x_dual, y_dual, info_dual] = sedumi(-model_dual.A, model_dual.b,model_dual.C,  model_dual.K);
%res_dual_unc = convert_sedumi2mosek(-model_dual.A',model_dual.b,model_dual.C,model_dual.K);
save('LR_120.mat', '-append', 'model_dual_cons_trans')
%save('LR_120.mat', '-append', 'res_dual_unc', 'model_dual_unc')