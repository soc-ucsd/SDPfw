N = 6;

[Mi, L] = tri_indexer(N);

cl = {[1,2],  [1, 4, 6], [3, 4, 5], [5, 6]};

s = cellfun(@length, cl);

i_curr = [];
j_curr = [];
v_curr = []; 

i_free = [];


Count_free = 0;
PSD_off = 100;
Count_rel = 1;
Count_psd = 0;
for k  = 1:length(cl)
    cl_curr = cl{k};
    N_curr = length(cl{k});
    Mi_curr = Mi(cl_curr, cl_curr);
    i_free = [i_free; Count_free + tri_vector(Mi_curr)];
    
    
    si_curr = reshape(Count_psd + (1:N_curr^2),  N_curr, N_curr);
    %vi2 = 
    for i = 1:N_curr
        for j = i:N_curr
            if i==j
                i_curr = [i_curr; si_curr(i, j)];
                j_curr = [j_curr; Count_rel];
                v_curr = [v_curr; -1];                
            else
                i_curr = [i_curr; si_curr(i, j); si_curr(j, i )];
                j_curr = [j_curr; Count_rel; Count_rel];
                v_curr = [v_curr; -0.5; -0.5];                
            end 
            Count_rel = Count_rel + 1;
        end                
    end  
    Count_psd = Count_psd + N_curr^2;
end

Count_rel = Count_rel - 1;
A_free_rel =  sparse(i_free, 1:Count_rel, ones(Count_rel, 1))';

A_s_rel = sparse(i_curr, j_curr, v_curr)';

A = [A_free_rel A_s_rel];