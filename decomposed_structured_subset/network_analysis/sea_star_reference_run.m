head_size = [0, 1, (1:14)*5];

Hout_list = zeros(length(head_size), 1);
time_list = zeros(length(head_size), 1);

outname = 'reference_results.mat';

QUIET = 0;

for k = 1:length(head_size)
% for k = 1:1
    curr_head = head_size(k);
    
    curr_name = strcat('sea_star_Hinf0_wide_med_', num2str(curr_head), '\\sea_star.mat');
    varlist = {'model', 'model_dense', 'n'};
    load(curr_name, varlist{:});
  
    N_agents = length(n);
    N_states = sum(n);
    
    [Hout, time_solve, sdp_opt] = sea_star_reference(model, QUIET);
    
    Hout_list(k) = Hout;
    time_list(k) = time_solve;
    
    save(outname, 'Hout_list', 'time_list', 'k', 'head_size')
    
    [k, Hout, time_solve];
    
    keyboard
    
end