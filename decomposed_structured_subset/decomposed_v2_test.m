%load('LR12_trans_box_sos.mat', 'model')

%model = model_cons_trans;
%model_new = struct;
model   = model_c;
%cones = cone_list(model.K.s, 10, "dd");
[model_new.A, model_new.b, model_new.C, model_new.K, model_new.info] =  ...
    decomposed_subset(model.A', model.b, model.c, model.K, cones);