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

[sdp_optP, cone_validP] = check_sdp_opt(xP ,yP_d, model.A, model.b, model.C, model.K, cones, 0);
prob = convert_sedumi2mosek(modelP.A, modelP.b, modelP.C, modelP.K);
[r, res] = mosekopt('minimize', prob);

[modelD.A, modelD.b, modelD.C, modelD.K, infoD] = ...
    decomposed_subset(model.A,model.b,model.C, model.K,cones, 1);

[xD_d, yD_d, info_solveD] = sedumi(modelD.A, modelD.b, modelD.C, modelD.K);
xD = decomposed_recover(xD_d, infoD);

XD1 = reshape(xD(2 + (1:9)), 3, 3);
XD2 = reshape(xD( 2 + 9 + (1:4)), 2, 2);


%When I do this to center X in DD* and get an initial feasible point for a 
%change of basis scheme, the answer I get is has X in DD, and the
%corresponding dual variable Z = 0. This is greatly unexpected. If instead
%I was working with a dual barrier to center Z, then a change of basis
%scheme could occur. What is strange is that we usually use dual barriers
%to center X. I don't think the dual change of basis is feasible without
%dualizing the problem and adding new variables.
[sdp_optD, cone_validD] = check_sdp_opt(xD ,yD_d, model.A, model.b, model.C, model.K, cones, 1);


[xD_0, yD_0, info_solveD_0] = sedumi(modelD.A, modelD.b, zeros(size(model.C)), modelD.K);
xD0 = decomposed_recover(xD_0, infoD);

XD10 = reshape(xD0(2 + (1:9)), 3, 3);
XD20 = reshape(xD0( 2 + 9 + (1:4)), 2, 2);

%z = recover_opt_dual


