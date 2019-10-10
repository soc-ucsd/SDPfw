%testing to make sure that the decomposed structured subset scheme works
rng(500, 'twister')

m = 100;
nBlk = 8;
BlkSize = 10;
ArrowHead = 10;


[model, model_split] = blockArrowSplit(m,nBlk,BlkSize,ArrowHead);


SP = spones(spones(model.c) + sparse(sum(spones(model.At),2)));  % vector of 1s and 0s
mask = reshape(SP, model.K.s, model.K.s);
%spy(mask)


%% SDP optimization


%now try the chordal decomposition
parCoLO.domain    = 1;  % dConvCliqueTree  ---> equalities 
parCoLO.range     = 2;   % rConvMatDecomp   ---> equalities 
parCoLO.EQorLMI   = 1; % CoLOtoEQform     ---> equality standard form
parCoLO.SDPsolver = []; % CoLOtoEQform     ---> equality standard form       
parCoLO.quiet     = 1; % Some peace and quiet       
J.f = length(model.b);

[~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(model.At',model.b,model.c,model.K,J,parCoLO); 

[x_psd, ~, ~] = sedumi(model.At, model.b, model.c, model.K);
[x_cpsd_G, ~, ~] = sedumi(LOP.A', LOP.b, LOP.c, LOP.K);
[x_cpsd, ~, ~] = sedumi(model_split.At, model_split.b, model_split.c, model_split.K);


%% Block Factor Width optimization

%Partition = 2;
Partition = Inf;
opts.bfw  = 1;
opts.nop  = Partition;
opts.socp = 1;   % second-order cone constraints

%SDD
[A_sdd, b_sdd, c_sdd, K_sdd, ~] = decomposed_subset(model.At',model.b,model.c,model.K,'sdd');

[x_sdd, ~, ~] = sedumi(A_sdd, b_sdd, c_sdd, K_sdd);

%SDD_clique
[A_csdd_G, b_csdd_G, c_csdd_G, K_csdd_G, ~] = decomposed_subset(LOP.A,LOP.b,LOP.c,LOP.K,'sdd');
[x_csdd_G, ~, ~] = sedumi(A_csdd_G, b_csdd_G, c_csdd_G, K_csdd_G);

[A_csdd, b_csdd, c_csdd, K_csdd, ~] = decomposed_subset(model_split.At,model_split.b,model_split.c,model_split.K,'sdd');
[x_csdd, ~, ~] = sedumi(A_csdd, b_csdd, c_csdd, K_csdd);

%DD
[A_dd, b_dd, c_dd, K_dd, ~] = decomposed_subset(model.At',model.b,model.c,model.K,'dd');

[x_dd, ~, ~] = sedumi(A_dd, b_dd, c_dd, K_dd);

%SDD_clique
[A_cdd_G, b_cdd_G, c_cdd_G, K_cdd_G, ~] = decomposed_subset(LOP.A,LOP.b,LOP.c,LOP.K,'dd');
[x_cdd_G, ~, ~] = sedumi(A_cdd_G, b_cdd_G, c_cdd_G, K_cdd_G);

[A_cdd, b_cdd, c_cdd, K_cdd, ~] = decomposed_subset(model_split.At,model_split.b,model_split.c,model_split.K,'dd');
[x_cdd, ~, ~] = sedumi(A_cdd, b_cdd, c_cdd, K_cdd);

%variable extraction is going to be more difficult

cost_sdd  = c_sdd'*x_sdd;
cost_csdd = c_csdd'*x_csdd;
cost_csdd_G = c_csdd_G'*x_csdd_G;

cost_dd  = c_dd'*x_dd;
cost_cdd = c_cdd'*x_cdd;
cost_cdd_G = c_cdd_G'*x_cdd_G;

cost_psd  = model.c'*x_psd;
cost_cpsd = model_split.c'*x_cpsd;
cost_cpsd_G = LOP.c'*x_cpsd_G;

Cost_Matrix = [cost_dd      cost_cdd_G cost_cdd;
               cost_sdd     cost_csdd_G cost_csdd;
               cost_psd     cost_cpsd_G cost_cpsd];
%dcost_sdd = cost_sdd_cl - cost_sdd_orig;

save('block_arrow.mat', 'model', 'model_split', 'LOP', 'cliqueDomain', 'Cost_Matrix')