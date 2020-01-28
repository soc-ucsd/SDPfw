%load('LR12_trans_box_sos.mat', 'model')
load('LR120_box_1_2.mat', 'model_cons_trans')
model = model_cons_trans;
%model_new = struct;
%model   = model_c;
%cones = cone_list(model.K.s, 12, dd);
cones = cone_list(model.K.s, 0, 2);
[model_new.A, model_new.b, model_new.C, model_new.K, model_new.info] =  ...
    decomposed_subset(model.A', model.b, model.c, model.K, cones);

 prob = convert_sedumi2mosek(model_new.A, model_new.b, model_new.C, model_new.K);
 
% [r, res] = mosekopt('minimize', prob)