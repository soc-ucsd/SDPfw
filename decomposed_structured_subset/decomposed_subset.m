function [Anew, bnew, cnew, Knew, info] = decomposed_subset(A,b,c,K,cones, dual)
%DECOMPOSED_SUBSET Reformulate a decomposed SDP into structured subsets
%   Each K.s clique is replaced by a cones in the list 'cones'
% Input data
%       A, b, c, K are SDP data in seudmi form
%       cones is a cell full of cones, one per entry in K.s.
%           like {'dd', 'dd', 'psd'}. If it is a string instead of a cell,
%           then all cliques get the same cones.
% Output data 
%       Anew, bnew, cnew, Knew, new SDP data in sedumi form
%       info.recover    A function that recovers output x into the desired
%           vector, immersion of structured subset into PSD.


%% Input check
 
    if nargin < 6
        dual = 0;
    end

    if size(A,1) ~= length(b) 
        A = A';
    end
    if ~isfield(K,'f') || isempty(K.f) 
        K.f = 0;
    end
    if ~isfield(K,'l') || isempty(K.l) 
        K.l = 0;
    end
    if ~isfield(K,'q') || isempty(K.q) 
        K.q = [];
    end

    m = length(b);
    
    %% Non PSD part
    %A_free_lin = A(:,1:K.f+K.l);
    A_free     = A(:, 1:K.f);
    A_lin      = A(:, K.f + (1:K.l));
    A_quad     = A(:, K.f + K.l + (1:K.q)); 
    %c_free_lin = c(1:K.f+K.l, :);
    c_free     = c(1:K.f);
    c_lin      = c(K.f + (1:K.l));
    c_quad     = c(K.f + K.l + (1:K.q), :);
    
    numPSD = length(K.s);
    
    
    %new cones
    Knew.f = K.f;
    Knew.l = K.l;
    Knew.q = K.q;
    Knew.s = [];
    
    %set up cones array
    if ~iscell(cones)
        cone_str = cones;
        cones = cell(numPSD, 1);
        for i = 1:numPSD
            cones{i} = cone_str;
        end
    end
    
    %define new matrices
    
    Count = K.f + K.l + sum(K.q);
    Count_free= 0;
    Count_linf= 0;
    Count_lin = 0;
    Count_psd = 0;
    Count_psdf= 0;
    Count_prev= 0;
        
    inew_linf = [];
    jnew_linf = [];
    vnew_linf = []; 
    
    inew_lin = [];
    jnew_lin = [];
    vnew_lin = [];    
    
    
    iprev_psd = [];
    jprev_psd = [];
    vprev_psd = [];
    Kprev_s = [];
    
    inew_psdf = [];
    jnew_psdf = [];
    vnew_psdf = [];
    
    
    inew_psd = [];
    jnew_psd = [];
    vnew_psd = [];
    
    %A_rel_free = {};
    A_rel_linf = {};
    A_rel_lin  = {};
    A_rel_psdf = {};
    A_rel_psd  = {};
    
    cnew_linf = [];
    cnew_lin  = [];        
    cprev_psd = [];
    cnew_psdf = [];
    cnew_psd  = [];
    %offset in current x for where the semidefinite variables are
            
    info_dd = {};
    info_non_dd = {};

    for i = 1:numPSD
        K_temp = struct;        
        
        cone_curr = cones{i};
        
        if strcmp(cone_curr, 'sdd')
            cone_curr = 1;
        end
        
        info_curr = struct;
        
        Ksi = K.s(i);
         
        info_curr.ind_orig = Count + (1:Ksi^2);
        
        A_curr = A(:, Count + (1:Ksi^2));
        c_curr = c(Count + (1:Ksi^2), :);
%         if Ksi == 1
%             %nonnegative entry, self dual
%             %Anew_lin = [Anew_lin A_curr];
%             [i_curr, j_curr, v_curr] = find(A_curr);
%             
%             j_curr = j_curr + Count_lin;
%             
%             %figure out how to preallocate this
%             inew_lin = [inew_lin; i_curr];
%             jnew_lin = [jnew_lin; j_curr];
%             vnew_lin = [vnew_lin; v_curr];                        
%             
%             cnew_lin = [cnew_lin; c_curr];
%             
%             info_curr.num_var = 1;
%             is_dd = 0;
%             
%             Knew.l = Knew.l + 1;
%             Count_lin = Count_lin + 1;
        %elseif strcmp('dd', cone_curr)
        if strcmp('dd', cone_curr)
            %DD
            K_temp.s = Ksi;
            if dual
                %DD star
                [~, ~, ~, ~, info_dds] = ...
                    dd_star_convert(A_curr, b, c_curr, K_temp, 1);
                
                A_dd_free = info_dds.A_dd_free ;
                A_dd_lin  = info_dds.A_dd_lin;                
                
                Knew.f = Knew.f + Ksi*(Ksi+1)/2;
                
                %New free variables
                [i_curr, j_curr, v_curr] = find(A_dd_free);

                j_curr = j_curr + Count_linf;
                
                inew_linf = [inew_linf; i_curr];
                jnew_linf= [jnew_linf; j_curr];
                vnew_linf= [vnew_linf; v_curr]; 

                
                %Indexing?
                info_curr.ind_dual = K.f + Count_linf + reshape(tri_indexer(Ksi),[], 1);
                
                Count_linf = Count_linf + Ksi*(Ksi+1)/2;         
                %Count_free = Count_free + Ksi*(Ksi+1)/2;
                info_curr.num_var = Ksi*(Ksi+1)/2;
                
                
                [i_curr, j_curr, v_curr] = find(A_dd_lin);

                j_curr = j_curr + Count_lin;

                %New nonnegative variables
                inew_lin = [inew_lin; i_curr];
                jnew_lin = [jnew_lin; j_curr];
                vnew_lin = [vnew_lin; v_curr]; 

                Count_lin = Count_lin + Ksi^2;
                Knew.l = Knew.l + Ksi^2;
                
                
                %Additional relationships (dual inequalities)
                A_rel_linf = cat(1,A_rel_linf, info_dds.A_rel_free);
                A_rel_lin = cat(1,A_rel_lin, info_dds.A_rel_lin);
                cnew_linf = [cnew_linf; info_dds.c_dd_free];
                cnew_lin   = [cnew_lin;  info_dds.c_dd_lin];
        
                
            else
                [A_dd_curr, ~, c_dd_curr, K_dd_curr, info_curr] = ...
                    dd_convert(A_curr, b, c_curr, K_temp);

                %info_curr.ind = info_curr.ind{1} + Count;
                info_curr.rays = info_curr.rays{1};
                info_curr.num_var = K_dd_curr.l;            

                %Anew_lin = [Anew_lin A_dd_curr];


                [i_curr, j_curr, v_curr] = find(A_dd_curr);

                j_curr = j_curr + Count_lin;

                %figure out how to preallocate this
                inew_lin = [inew_lin; i_curr];
                jnew_lin = [jnew_lin; j_curr];
                vnew_lin = [vnew_lin; v_curr]; 

                Count_lin = Count_lin + Ksi^2;
                Knew.l = Knew.l + Ksi^2;

                cnew_lin = [cnew_lin; c_dd_curr];                 
            end
            
            is_dd = 1;       
        elseif isnumeric(cone_curr)  && (cone_curr > 0)           
            %SDD
            opts.bfw  = 1;
            %opts.socp = 1;   % second-order cones constraints
            opts.socp = 0;   % No SOCP constraints, they are bugged
            opts.block = cone_curr;
            
            
            
            if dual
                %Factor Width star
                opts.keep_split = 1;
                opts.dual = 1;
                K_temp.s = Ksi;
                [~, ~, ~, ~, info_fw] = ...
                    factorwidth(A_curr, b, c_curr, K_temp, opts);                            
                
                A_fw_free = info_fw.A_fw_free ;
                A_fw_psd  = info_fw.A_fw_psd;                
                
                Knew.f = Knew.f + Ksi*(Ksi+1)/2;
                
                %New free variables
                [i_curr, j_curr, v_curr] = find(A_fw_free);

                j_curr = j_curr + Count_psdf;
                
                inew_psdf = [inew_psdf; i_curr];
                jnew_psdf= [jnew_psdf; j_curr];
                vnew_psdf= [vnew_psdf; v_curr]; 

                
                %redo this
                info_curr.ind_dual = K.f + Count_psdf + reshape(tri_indexer(Ksi),[], 1);
                
                Count_psdf = Count_psdf + Ksi*(Ksi+1)/2;
                %Count_free = Count_free + Ksi*(Ksi+1)/2;         
                
                info_curr.num_var = Ksi*(Ksi+1)/2;
                
                num_s_new = length(info_fw.c_fw_psd);
                
                [i_curr, j_curr, v_curr] = find(A_fw_psd);

                j_curr = j_curr + Count_psd;

                %New nonnegative variables
                inew_psd = [inew_psd; i_curr];
                jnew_psd = [jnew_psd; j_curr];
                vnew_psd = [vnew_psd; v_curr]; 

                Count_psd = Count_psd + num_s_new;
                Knew.s = [Knew.s;  info_fw.new_cone.s];
                
                
                %Additional relationships (dual inequalities)
                A_rel_psdf = cat(1,A_rel_psdf, info_fw.A_rel_free);
                A_rel_psd  = cat(1,A_rel_psd,  info_fw.A_rel_psd);
                
                cnew_psdf  = [cnew_psdf; info_fw.c_fw_free];
                cnew_psd   = [cnew_psd;  info_fw.c_fw_psd];
        
                
            else
                K_temp.s = Ksi;
                %TODO: add indexing/recovery for sdd_info
                opts.keep_split = 0;
                opts.dual = 0;
                [A_curr, ~, c_curr, K_curr, info_curr] = ...
                    factorwidth(A_curr, b, c_curr, K_temp, opts);            

                [i_curr, j_curr, v_curr] = find(A_curr);

                j_curr = j_curr + Count_psd;

                %figure out how to preallocate this
                inew_psd = [inew_psd; i_curr];
                jnew_psd = [jnew_psd; j_curr];
                vnew_psd = [vnew_psd; v_curr]; 

                Count_psd = Count_psd + sum(K_curr.s.^2);

                
                cnew_psd = [cnew_psd; c_curr];


                Knew.s = [Knew.s; K_curr.s];
                info_curr.num_var = length(c_curr);
            end
            is_dd = 0;
        else
            %PSD, self dual
            %Anew_psd = [Anew_psd A_curr];
            
            [i_curr, j_curr, v_curr] = find(A_curr);
            
            j_curr = j_curr + Count_prev;
            
            %figure out how to preallocate this
            iprev_psd = [iprev_psd; i_curr];
            jprev_psd = [jprev_psd; j_curr];
            vprev_psd = [vprev_psd; v_curr]; 
            
            cprev_psd = [cprev_psd; c_curr];
            info_curr.num_var = length(c_curr);
            is_dd = 0;
            Count_prev = Count_prev+ Ksi^2;
            Kprev_s = [Kprev_s; Ksi];
        end
        info_curr.ind_orig = Count + (1:Ksi^2);
        
        if is_dd
            info_dd{end+1} = info_curr;
        else
            info_non_dd{end+1} = info_curr;
        end
        Count = Count + Ksi^2;
    end

    %output results
    
    Knew.s = [Kprev_s; Knew.s];
    %Anew_free = sparse(inew_free, jnew_free, vnew_free, m, Count_free);
    Anew_free_lin = sparse(inew_linf, jnew_linf, vnew_linf, m, Count_linf);
    Anew_free_psd = sparse(inew_psdf, jnew_psdf, vnew_psdf, m, Count_psdf);
    
    Anew_lin  = sparse(inew_lin,  jnew_lin,  vnew_lin, m, Count_lin);
    Aprev_psd = sparse(iprev_psd,  jprev_psd,  vprev_psd, m, sum(Kprev_s.^2));
    Anew_psd  = sparse(inew_psd,  jnew_psd,  vnew_psd, m, Count_psd);
    
    Anew = [A_free Anew_free_lin Anew_free_psd A_lin Anew_lin A_quad Aprev_psd Anew_psd];
    
    cnew = [c_free; cnew_linf; cnew_psdf; c_lin; cnew_lin; c_quad; cprev_psd; cnew_psd];
    
     
    Knew.q(Knew.q == 0) = [];
    
    
        %relations between data
    if dual
        
        for k = 1:length(info_non_dd)
            if isfield(info_non_dd{k}, 'ind_dual')
                info_non_dd{k}.ind_dual = info_non_dd{k}.ind_dual + Count_linf;
            end
        end
        
        %A_rel_free_diag = cell_blkdiag(A_rel_free);
        A_rel_free_lin_diag = cell_blkdiag(A_rel_linf);        
        A_rel_lin_diag  = cell_blkdiag(A_rel_lin);
        num_rel_lin =  size(A_rel_free_lin_diag, 1);
        
        
        A_rel_free_psd_diag = cell_blkdiag(A_rel_psdf);        
        A_rel_psd_diag  = cell_blkdiag(A_rel_psd);
        num_rel_psd =  size(A_rel_free_psd_diag, 1);
    
        %A_rel_lin_diag =  -speye(num_rel);
        
        A_rel_lin_only = [sparse(num_rel_lin, K.f), A_rel_free_lin_diag, sparse(num_rel_lin, length(cnew_psdf) + K.l), ...
            A_rel_lin_diag, sparse(num_rel_lin, sum(K.q) + sum(Knew.s.^2))];
            
        A_rel_psd_only = [sparse(num_rel_psd, K.f), sparse(num_rel_psd, length(cnew_linf)),A_rel_free_psd_diag, ...
            sparse(num_rel_psd, Knew.l + sum(K.q)), sparse(num_rel_psd, sum(Kprev_s.^2)),  A_rel_psd_diag];
            
        %cnew = [c_free; cnew_linf; cnew_psdf; c_lin; cnew_lin; c_quad; cprev_psd; cnew_psd];    
        Anew = [Anew; A_rel_lin_only; A_rel_psd_only];    
        
        %This will not work with mixes of  psd and FW*. Handle this 
        %A_rel = [sparse(num_rel, K.f), A_rel_free_diag, sparse(num_rel, K.l), ...
        %    A_rel_lin_diag, sparse(num_rel, sum(K.q) + sum(Kprev_s.^2)), A_rel_psd_diag];
        
        %Anew = [Anew; A_rel];
        num_rel = num_rel_lin + num_rel_psd;
        bnew = [b; sparse(num_rel, 1)] ;
    else
        bnew = b;
        num_rel = 0;
    end
    
    
    %info = [];
    info.dual = dual;
    info.dd = info_dd;    
    info.num_rel = num_rel;
    info.count_free = Count_free;
    info.non_dd = info_non_dd;
    info.prev_cone = K;
    info.count_non_dd = Knew.f + Knew.l + sum(Knew.q);
    info.prev_var = length(c);
    info.dual = dual;
    
end

