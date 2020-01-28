function [Anew, bnew, cnew, Knew, info] = decomposed_subset(A,b,c,K,cones)
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

    %% Non PSD part
    A_free_lin = A(:,1:K.f+K.l);
    A_quad     = A(:, K.f + K.l + (1:K.q)); 
    c_free_lin = c(1:K.f+K.l, :);
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
    
    
    %first pass through the PSD cones for allocation
    %offsets for the cones in the system
    %reformulated problem's linear variables will start at Count_l (for
    %example) 
    %Count_lin = K.f + K.l; 
    %Count_psd = K.f + K.l + sum(K.q);
    
%     new_l = 0;
%     new_s = [];
%     for i = 1:numPSD
%         Ksi = K.s(i);
%         if strcmp('sdd', cones{i}) 
%             cones{i} = 1;
%         end
%         
%         cone_curr = cones{i};
%         
%         if Ksi == 1
%             %singleton PSD  block ==> linear constraint
%             Knew.l = Knew.l + 1;
%             %Count_psd = Count_psd + 1;
%         elseif strcmp('dd', cone_curr)
%             %diagonally dominant
%             %Count_psd = Count_psd + Ksi^2;
%             Knew.l = Knew.l + Ksi^2;
%         elseif isnumeric(cone_curr)  && (cone_curr > 0)    
%             %block-factor width
%             
%             %Copy and paste from factorwidth.m
%             %I know, DRY and all. still.
%             %number of partition
%             nop = max(floor(Ksi/cone_curr), 1);
%             
%             SizeU = ceil(K.s(PSDind)/nop);
%             SizeL = floor(K.s(PSDind)/nop);
%             if SizeU == SizeL
%                 alpha = ones(nop,1)*SizeU;
%             else
%                 x = (K.s(PSDind) - SizeL*nop)./(SizeU-SizeL);
%                 alpha = [ones(x,1)*SizeU;ones(nop-x,1)*SizeL];
%             end
%             Knew.s = [Knew.s; alpha];
%         else
%             %PSD
%             Knew.s = [Knew.s; Ksi];
%         end                
%     end
    
    
    %Nvars = Knew.f + Knew.l + Knew.q + sum(Knew.s.^2);
    
    %Anew_lin  = [];
    %Anew_quad = [];
    %Anew_psd  = [];
    
    %define new matrices
    
    Count = K.f + K.l + sum(K.q);
    
    Count_lin = 0;
    Count_psd = 0;
    
    inew_lin = [];
    jnew_lin = [];
    vnew_lin = [];    
    
    inew_psd = [];
    jnew_psd = [];
    vnew_psd = [];
    
    cnew_lin  = [];
    %cnew_quad = [];
    cnew_psd  = [];
    %offset in current x for where the semidefinite variables are
            
    info_dd = {};
    info_non_dd = {};

    for i = 1:numPSD
        K_temp = struct;        
        
        cone_curr = cones{i};
        info_curr = struct;
        
        Ksi = K.s(i);
         
        info_curr.ind_orig = Count + (1:Ksi^2);
        
        A_curr = A(:, Count + (1:Ksi^2));
        c_curr = c(Count + (1:Ksi^2), :);
        if Ksi == 1
            %nonnegative entry, self dual
            %Anew_lin = [Anew_lin A_curr];
            [i_curr, j_curr, v_curr] = find(A_curr);
            
            j_curr = j_curr + Count_lin;
            
            %figure out how to preallocate this
            inew_lin = [inew_lin; i_curr];
            jnew_lin = [jnew_lin; j_curr];
            vnew_lin = [vnew_lin; v_curr]; 
            
            Count_lin = Count_lin + Ksi^2;
            
            cnew_lin = [cnew_lin; c_curr];
            
            info_curr.num_var = 1;
            is_dd = 0;
            
            Knew.l = Knew.l + 1;
            Count_lin = Count_lin + 1;
        elseif strcmp('dd', cone_curr)
            %DD
            K_temp.s = Ksi;
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
            
            is_dd = 1;            
        elseif isnumeric(cone_curr)  && (cone_curr > 0)           
            %SDD
            opts.bfw  = 1;
            %opts.socp = 1;   % second-order cones constraints
            opts.socp = 0;   % No SOCP constraints, they are bugged
            opts.block = cone_curr;

            
            K_temp.s = Ksi;
            %TODO: add indexing/recovery for sdd_info
            [A_curr, ~, c_curr, K_curr, info_curr] = ...
                factorwidth(A_curr, b, c_curr, K_temp, opts);            

            
%             if opts.socp && (cone_curr == 1)
%                 K_curr.q(K_curr.q == 0) = [];
%                 Knew.q = [Knew.q K_curr.q];
%                 
%                 Anew_quad = [Anew_quad A_curr];
%                 cnew_quad = [cnew_quad; c_curr];
%                 
%             else
            [i_curr, j_curr, v_curr] = find(A_curr);
            
            j_curr = j_curr + Count_psd;
            
            %figure out how to preallocate this
            inew_psd = [inew_psd; i_curr];
            jnew_psd = [jnew_psd; j_curr];
            vnew_psd = [vnew_psd; v_curr]; 
            
            Count_psd = Count_psd + sum(K_curr.s.^2);

            %Anew_psd = [Anew_psd K_curr];
            cnew_psd = [cnew_psd; c_curr];
            
            
            Knew.s = [Knew.s; K_curr.s];
            info_curr.num_var = length(c_curr);
            is_dd = 0;
        else
            %PSD, self dual
            %Anew_psd = [Anew_psd A_curr];
            
            [i_curr, j_curr, v_curr] = find(A_curr);
            
            j_curr = j_curr + Count_psd;
            
            %figure out how to preallocate this
            inew_psd = [inew_psd; i_curr];
            jnew_psd = [jnew_psd; j_curr];
            vnew_psd = [vnew_psd; v_curr]; 
            
            cnew_psd = [cnew_psd; c_curr];
            info_curr.num_var = length(c_curr);
            is_dd = 0;
            Count_psd = Count_psd + Ksi^2;
            Knew.s = [Knew.s; Ksi];
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
    %Anew = [A_free_lin Anew_lin A_quad Anew_quad Anew_psd];
    
    Anew_lin = sparse(inew_lin, jnew_lin, vnew_lin, length(b) ,Count_lin);
    Anew_psd = sparse(inew_psd, jnew_psd, vnew_psd, length(b), Count_psd);
    
    Anew = [A_free_lin Anew_lin A_quad Anew_psd];
    
    cnew = [c_free_lin; cnew_lin; c_quad; cnew_psd];
    bnew = b;
    
    Knew.q(Knew.q == 0) = [];
    
    %info = [];
    
    info.dd = info_dd;
    info.count_dd = K.f+K.l;
    info.non_dd = info_non_dd;
    info.count_non_dd = Count_lin + sum(K.q);
    info.prev_var = length(c);
    
end

