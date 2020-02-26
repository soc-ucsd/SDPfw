load('sos_quartic_200.mat', 'model_split');
%load('sos_quartic_200.mat', 'model');
model = model_split;
modelD = struct;
model.c =  model.C;
model = rmfield(model, 'c');

%  [modelD.A, modelD.b, modelD.C, modelD.K, info] = ...
%      dd_star_convert(model.A,model.b,model.C, model.K);
% 
%  

opt.bfw = 1;
opt.block = 1 ;
opt.dual = 1;
% 
% [modelF.A, modelF.b, modelF.C, modelF.K, infoF] = ...
%     factorwidth(model.A,model.b,model.C, model.K, opt);
 
[modelS.A, modelS.b, modelS.C, modelS.K, infoS] = ...
    decomposed_subset(model.A,model.b,model.C, model.K, {'dd','psd'}, 0);


%[xf, yf, info_solve] = sedumi(modelF.A, modelF.b, modelF.C, modelF.K);
%  
% pars.fid = 0;
% % [xd, yd, info_D] = sedumi(modelD.A, modelD.b, modelD.C, modelD.K, pars);
[xs, ys, info_solve] = sedumi(modelS.A, modelS.b, modelS.C, modelS.K);
x = decomposed_recover(xs, infoS);

X1 = reshape(x(2 + (1:9)), 3, 3);
X2 = reshape(x( 2 + 9 + (1:4)), 2, 2);
