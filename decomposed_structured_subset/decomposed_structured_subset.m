%testing to make sure that the decomposed structured subset scheme works
rng(500, 'twister')

m = 40;
nBlk = 10;
BlkSize = 6;
ArrowHead = 4;
r = 6;

[At,b,c,K, ~] = blockArrowRank(m,nBlk,BlkSize,ArrowHead, r);


SP = spones(spones(c) + sparse(sum(spones(At),2)));  % vector of 1s and 0s
mask = reshape(SP, K.s, K.s);
%spy(mask)


%% SDP optimization
[x, ~, ~] = sedumi(At, b, c, K);

%now try the chordal decomposition
parCoLO.domain    = 1;  % dConvCliqueTree  ---> equalities 
parCoLO.range     = 2;   % rConvMatDecomp   ---> equalities 
parCoLO.EQorLMI   = 1; % CoLOtoEQform     ---> equality standard form
parCoLO.SDPsolver = []; % CoLOtoEQform     ---> equality standard form       
parCoLO.quiet     = 1; % Some peace and quiet       
J.f = length(b);

[~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(At',b,c,K,J,parCoLO); 

[x_cl, ~, ~] = sedumi(LOP.A', LOP.b, LOP.c, LOP.K);

cost_orig = c'*x;
cost_cl = LOP.c'*x_cl;

dcost = cost_cl - cost_orig;

%% Block Factor Width optimization

%Partition = 2;
Partition = Inf;
opts.bfw  = 1;
opts.nop  = Partition;
opts.socp = 1;   % second-order cone constraints

%standard block factor width
[A_fw, b_fw, c_fw, K_fw, info_fw_decomp] = factorwidth(At',b,c,K,opts);

[x_fw, ~, ~] = sedumi(A_fw, b_fw, c_fw, K_fw);

%block factor width over cliques
[A_cl_fw, b_cl_fw, c_cl_fw, K_cl_fw, info_cl_fw_decomp] = factorwidth(LOP.A,LOP.b,LOP.c,LOP.K,opts);

[x_cl_fw, ~, ~] = sedumi(A_cl_fw, b_cl_fw, c_cl_fw, K_cl_fw);

%variable extraction is going to be more difficult

cost_fw_orig = c_fw'*x_fw;
cost_fw_cl = c_cl_fw'*x_cl_fw;
dcost_fw = cost_fw_cl - cost_fw_orig;