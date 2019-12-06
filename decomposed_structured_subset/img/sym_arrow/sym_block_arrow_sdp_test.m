rng(10, 'twister')
visualize = 1;
%number of matrices
m =  40;
ArrowHead = 10;
BlkSize = 10;

% nBlk = 6;
% orbits = {[1,2,3,4],5, 6};
% 
% visualize = 1;
% 

%nBlk = 8;
%orbits = {[1,2,3,4], 5, 6, [7,8]};
nBlk = 4;
%orbits =  {1,2,3,4};
orbits =  {[1,2],3,4};

orbit_ind = cellfun(@(x) (x-1)*BlkSize +  (1:BlkSize)', orbits , 'UniformOutput', false);
arrow_ind = nBlk*BlkSize + (1:ArrowHead);

N = nBlk*BlkSize + ArrowHead;

% nBlk = 10;
% orbits = {[1,2,3], 4, [5, 6], [7,  8,9, 10]};
%orbits = {[1,3,4,9], [2,10], [5,8],  6, 7};
%orbits = {1:10};
[model, model_sym] = sym_block_arrow_sdp(m, nBlk, BlkSize, ArrowHead, orbits);

[x,y,info] = sedumi(model.At, model.b, model.c, model.K, model.pars);
[xs,ys,info] = sedumi(model_sym.At, model_sym.b, model_sym.c, model_sym.K, model_sym.pars);

cost = model.c'*x;
cost_sym = model_sym.c'*xs;