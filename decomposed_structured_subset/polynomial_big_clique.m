load('polynomial_big_clique.mat')
% 
% N = 64;
% b = N/4;
% 
% 
 cones = {'dd', 'sdd', 2, 10, 30, 'psd'}; 
 lc = length(cones);
ls = length(model_csp.K.s);
cone_list = cell(ls, 1);

for k = 1:ls
    cone_list{k} = 'psd';
end

model_list = cell(lc, 2);

for i = 1:lc
    model_curr = struct;
    %all in same clique 
    [model_curr.A, model_curr.b, model_curr.C, model_curr.K, ~] = decomposed_subset(model_csp.A,model_csp.b,model_csp.C,model_csp.K, cones{i});
    model_curr.pars = model_csp.pars;
    model_list{i, 1} = model_curr;
    if sum(model_curr.K.q == 0) 
        keyboard
    end
    
    %last clique is different
    cone_list{end} = cones{i};
    [model_curr.A, model_curr.b, model_curr.C, model_curr.K, ~] = decomposed_subset(model_csp.A,model_csp.b,model_csp.C,model_csp.K, cone_list);
    %model_curr.K
    model_curr.pars = model_csp. pars;
    model_list{i, 2} = model_curr;  
    if sum(model_curr.K.q == 0) 
        keyboard
    end
        
end

save('polynomial_big_clique.mat', 'model_list', 'model_csp', 'cones', 'b', 'N')

RUN_DATA = 1;
%lc = length(cones);

if RUN_DATA
    cost_list = zeros(lc, 2);
    info_list = cell(lc, 2);
end

if RUN_DATA
    for i = 1:lc
        {i, 'same'}  
        mc = model_list{i,1};
        [x, y, info] = sedumi(mc.A, mc.b, mc.C, mc.K, mc.pars);
        cost_list(i, 1) = mc.C'*x;
        info_list{i, 1} = info;
        
        {i, 'diff'}  
        mc = model_list{i,2};
        [x, y, info] = sedumi(mc.A, mc.b, mc.C, mc.K, mc.pars);
        cost_list(i, 2) = mc.C'*x;                
        info_list{i, 2} = info;
    end
end