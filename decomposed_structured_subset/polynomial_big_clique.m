PROCESS_DATA = 1;
RUN_DATA = 1;


if PROCESS_DATA
    %load('LR_120_sos.mat')
    load('LR_120_sos.mat')
    model_csp.pars.fid = 0;
    
    %cones = {'dd', 'sdd', 2, ceil(maxKs/8), ceil(maxKs/4), 'psd'}; 
    %cones = {'dd', 'sdd', 2, 10, 19, 'psd'};
    %cones = {'dd', 'sdd', 2, 3, 5, 10, 15, 30, 'psd'};
    cones = {'dd', 'sdd', 2, 3, 5, 11, 21, 40, 'psd'};
    %cones = {4, 7, 9};
    %cones = {6};
    %cones = {'dd', 5, 19};
    
    %cones = {'dd', 'sdd', 2, 3, 5,  10, 19, 'psd'}; 
    %cones = {'dd', 'sdd', 9, 'psd'};
    %cones = {'dd', 1, 2, 3, 5, 'psd'}; 
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
        model_curr.pars = model_csp.pars;
        model_list{i, 2} = model_curr;  
        if sum(model_curr.K.q == 0) 
            keyboard
        end

    end

    save('LR_120_sos_process.mat', 'model_list', 'model_csp', 'cones', 'b', 'N')

else
    load('LR_120_sos_process.mat')
end

lc = length(cones);
maxKs = max(model_csp.K.s);

if RUN_DATA
    cost_list = zeros(lc, 2);
    info_list = cell(lc, 2);
    time_list = zeros(lc, 2);
end

if RUN_DATA
    i = 4;
    for i = 1:lc
        {cones{i}, 'same'}  
        [time_list(i, 1), cost_list(i, 1), info_list{i, 1}, x1] = run_model(model_list{i, 1});
        
        {cones{i}, 'diff'}          
        [time_list(i, 2), cost_list(i, 2), info_list{i, 2}, x2] = run_model(model_list{i, 2});
        save('LR_120_sos_jan23.mat', 'time_list', 'cost_list', 'info_list', 'output_table')
    end
end
output_table = latex(vpa(sym([cost_list time_list]),2));


save('LR_120_sos_120_jan23.mat', 'time_list', 'cost_list', 'info_list', 'output_table')

%save('big_clique_results_fR2_120.mat', 'time_list', 'cost_list', 'info_list', 'output_table')
%save('big_clique_results_32.mat_0', 'time_list', 'cost_list', 'info_list')
