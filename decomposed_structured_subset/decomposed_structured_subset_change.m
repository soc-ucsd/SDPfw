%testing to make sure that the decomposed structured subset scheme works
rng(500, 'twister')
% 
SOLVE = 0;
PLOT = 1;

m = 80;
nBlk = 15;
BlkSize = 10;
ArrowHead = 10;

% m = 80;
% nBlk = 5;
% BlkSize = 10;
% ArrowHead = 5;
if SOLVE

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

LOP.K.f = 0;
LOP.K.l = 0;
LOP.K.q = [];


% mask_CoLO = sparse(size( mask, 1), size( mask, 2));
% for i = 1:cliqueDomain{1,1}.NoC
%     cli = cliqueDomain{1,1}.Set{i};
%     mask_CoLO(cli, cli) = 1;
% end 


% figure(1)
% hold on
% spy(mask_CoLO , 'm');
% spy(mask);
% hold off
% title('Block Arrow Sparsity + Fill-in', 'fontsize', 14)



LOP.At = LOP.A';


iter_max = 20;
eps_change = 1e-1;

%cone_list = {'sdd', 5, 10, 'psd'};
%cone_list = {'dd', 'sdd', 2, 5, 10, 'psd'};
cone = 5;


model_list = {model, LOP, model_split};

model_name = {'K', 'K(E_colo, ?)', 'K(E_sparse, ?)'};
Cost_Matrix = zeros(iter_max, length(model_list));
Time_Matrix = zeros(iter_max, length(model_list));
Info_Matrix = cell(iter_max, length(model_list));
X_Matrix = cell(iter_max+1, length(model_list));

x_fake = cell(length(model_list), 1);

Basis_Matrix = cell(iter_max+1, length(model_list));

%start with initial feasible point
for j = 1:length(model_list)
    model_curr = model_list{j};
    model_feas = model_curr;
    model_feas.c = zeros(size(model_feas.c));
    
    [cost, x_feas, info, time] = run_model(model_feas, 'psd');
    
    [model_feas_change, x_fake{j}] = basis_change(x_feas, model_curr);
    Basis_Matrix{1, j} = model_feas_change; 
    X_Matrix{1, j} = x_feas;
    %Basis_Matrix{1, i} = model_feas_change;
end
% Basis_Matrix{1, 1} = model;
% Basis_Matrix{1, 2} = LOP;
% Basis_Matrix{1, 3} = model_split;

%Full SDP
[cost0, x0, info0, time0] = run_model(LOP, 'psd');

%change of basis
for i = 1:iter_max 
    for j = 1:length(model_list)            
        {model_name{j}, i}
        curr_model = Basis_Matrix{i, j};
        
        [cost, x, info, time] = run_model(curr_model, cone);
        Cost_Matrix(i, j) = cost;        
        Info_Matrix{i, j} = info;
        Time_Matrix(i, j) = time;
        
        X_Matrix{i+1, j}    = x;
        %x_change = (1-eps_change)*x + eps_change*X_Matrix{1, j};
        x_change = (1-eps_change)*x + eps_change*x_fake{j};
        Basis_Matrix{i+1, j} = basis_change(x_change, Basis_Matrix{i,j});
    end 
end

save('block_arrow_change.mat', 'cone', 'cliqueDomain', ...
    'Cost_Matrix', 'Info_Matrix', 'Basis_Matrix', 'X_Matrix', 'cost0', 'x0', 'model_list')

end
if PLOT
figure(1)
clf
C = linspecer(5);
hold on
plot(1:iter_max, Cost_Matrix(:, 1), 's-', 'color', C(1, :), 'linewidth', 2)
plot(1:iter_max, Cost_Matrix(:, 2), 's-','color', C(2, :), 'linewidth', 2)
plot(1:iter_max, Cost_Matrix(:, 3), 's-', 'color', C(3, :), 'linewidth', 2)

xlabel('Iteration', 'FontSize', 14)
ylabel('Objective', 'FontSize', 14)

plot([1, iter_max], [cost0, cost0], 'k--', 'linewidth', 2)
hold off
legend_list = {'$B_5$', '$B_5(E_F, ?)$', '$B_5(E, ?)$', '$S_+$'};
legend(legend_list, 'Interpreter', 'latex', 'FontSize', 14)
title('Decomposed Change of Basis Comparision', 'FontSize', 18)

end
function  [cost, x, info, time] = run_model(model, cone)
    tic
    [A, b, c, K, info_new] = decomposed_subset(model.At,model.b,model.c,model.K, cone);
    time = toc;
    pars.fid = 0;
    [x, ~, info] = sedumi(A, b, c, K, pars);
    cost = c'*x;
    x = decomposed_recover(x, info_new);    
end
