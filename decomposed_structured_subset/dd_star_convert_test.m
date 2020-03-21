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
 
cones =  {'sdd','dd'};

[modelP.A, modelP.b, modelP.C, modelP.K, infoP] = ...
    decomposed_subset(model.A,model.b,model.C, model.K,cones, 0);

[xP_d, yP_d, info_solveP] = sedumi(modelP.A, modelP.b, modelP.C, modelP.K);


xP = decomposed_recover(xP_d, infoP);

XP1 = reshape(xP(2 + (1:9)), 3, 3);
XP2 = reshape(xP( 2 + 9 + (1:4)), 2, 2);

[sdp_optP, cone_validP] = check_sdp_opt(xP_d ,yP_d, model.A, model.b, model.C, model.K, cones, 0);
prob = convert_sedumi2mosek(modelP.A, modelP.b, modelP.C, modelP.K);
[r, res] = mosekopt('minimize', prob);
% [modelD.A, modelD.b, modelD.C, modelD.K, infoD] = ...
%     decomposed_subset(model.A,model.b,model.C, model.K,cones, 1);
% 
% [xD_d, yD_d, info_solveD] = sedumi(modelD.A, modelD.b, modelD.C, modelD.K);
% xD = decomposed_recover(xD_d, infoD);
% 
% XD1 = reshape(xD(2 + (1:9)), 3, 3);
% XD2 = reshape(xD( 2 + 9 + (1:4)), 2, 2);
% 
% [sdp_optD, cone_validD] = check_sdp_opt(xD ,yD_d, model.A, model.b, model.C, model.K, cones, 1);
%z = recover_opt_dual


