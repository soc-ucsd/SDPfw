rng(10, 'twister')
visualize = 1;
%number of matrices
m =  200;
ArrowHead = 20;
BlkSize = 20;

% nBlk = 6;
% orbits = {[1,2,3,4],5, 6};
% splits = [1 1 0];
% visualize = 1;


nBlk = 8;
orbits = {[1,2,3,4], 5, 6, [7,8]};
splits = [1, 0, 1, 0 ];


%nBlk = 4;
%orbits =  {1,2,3,4};
%orbits =  {[1,2],3,4};

orbit_ind = cellfun(@(x) (x-1)*BlkSize +  (1:BlkSize)', orbits , 'UniformOutput', false);
arrow_ind = nBlk*BlkSize + (1:ArrowHead);

N = nBlk*BlkSize + ArrowHead;

% nBlk = 10;
% orbits = {[1,2,3], 4, [5, 6], [7,  8,9, 10]};
%orbits = {[1,3,4,9], [2,10], [5,8],  6, 7};
%orbits = {1:10};
[model, model_sym] = sym_block_arrow_sdp(m, nBlk, BlkSize, ArrowHead, orbits, splits);

cone = 1;


% [x,y,info] = sedumi(model.At, model.b, model.c, model.K, model.pars);
% [xs,ys,info] = sedumi(model_sym.At, model_sym.b, model_sym.c, model_sym.K, model_sym.pars);
% 
% cost = model.c'*x;
% cost_sym = model_sym.c'*xs;

%psd cone
%[cost, ~, ~] = run_model (model, 'psd');
%[cost_sym, ~, ~] = run_model (model_sym, 'psd');


%sparsecolo

%now try the chordal decomposition
parCoLO.domain    = 1;  % dConvCliqueTree  ---> equalities 
parCoLO.range     = 2;   % rConvMatDecomp   ---> equalities 
parCoLO.EQorLMI   = 1; % CoLOtoEQform     ---> equality standard form
parCoLO.SDPsolver = []; % CoLOtoEQform     ---> equality standard form       
parCoLO.quiet     = 1; % Some peace and quiet       
%J.f = length(model.b);
J = [];

[~,~,~,cliqueDomain,cliqueRange,model_sp] = sparseCoLO(model.At',model.b,model.c,model.K,J,parCoLO); 
[~,~,~,cliqueDomain_sym,cliqueRange_sym,model_sym_sp] = sparseCoLO(model_sym.At',model_sym.b,model_sym.c,model_sym.K,J,parCoLO); 
model_sp.At = model_sp.A';
model_sym_sp.At = model_sym_sp.A';
[cost_sp, ~, ~] = run_model (model, 'psd');
[cost_sym_sp, ~, ~] = run_model (model_sym, 'psd');


%structured subsets
cone = 1;
[costK, infoK, timeK]        = run_model(model, cone);
[costK_sym, infoK_sym, timeK_sym]    = run_model(model_sym, cone);
[costK_sp,  infoK_sp , timeK_sp]     = run_model(model_sp, cone);
[costK_sym_sp, infoK_sym_sp, timeK_sym_sp] = run_model(model_sym_sp, cone);


%cost0_mat = [cost cost_sym; cost_sp costt_sym_sp]; 
cost_mat = [costK costK_sym; costK_sp costK_sym_sp]; 
time_mat = {timeK timeK_sym; timeK_sp timeK_sym_sp};
time_pose_mat = cellfun(@(x) x(1) , time_mat);
time_solve_mat = cellfun(@(x) x(2), time_mat);
time_total_mat = cellfun(@(x) x(1)+x(2), time_mat);
function  [cost, info, time] = run_model(model, cone)
    tic
    [A, b, c, K, ~] = decomposed_subset(model.At,model.b,model.c,model.K, cone);
    time_pose = toc;
    pars.fid = 0;
    tic
    [x, ~, info] = sedumi(A, b, c, K, pars);
    time_solve = toc;
    cost = c'*x;
    time  = [time_pose ,time_solve];
end