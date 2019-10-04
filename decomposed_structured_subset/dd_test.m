load('quartic_process_3m.mat')


model_dd = struct;
[model_dd.At, model_dd.b, model_dd.c, model_dd.K, model_dd.info] = ...
    dd_convert(model.At,model.b,model.c,model.K);

[x, y, info] = sedumi(model_dd.At, model_dd.b, model_dd.c, model_dd.K);

%Cmodel_dd = struct;