function [Anew, bnew, cnew, Knew, info] = factorwidth(A,b,c,K,opts)
%  Reformulating a primal SDP with a block factorwidth two cone
%
%       min_{x} c'*x
%               Ax = b
%               x \in K
%
% K can have K.f, K.l, K.q, K.s; 
%       Only replacing K.s with a block factor-width-two cone
%       and reformulating it into a standard SDP in the SeDuMi form
% Input data
%       A, b, c, K are SDP data in seudmi form
%       opts.bfw     1 or 0,  block factor-width-two decomposition
%       opts.nop     integer, number of blocks in the partion alpha
%       opts.size    alternative to nop, number of entries in each block
%       opts.socp    1 or 0,  reformualte 2 by 2 PSD cone with a second-order cone
%       opts.dual    1 or 0, whether this should be dual or primal block
%                    factorwidth two cone
% Output data 
%       Anew, bnew, cnew, Knew, new SDP data in sedumi form
%       info.Ech    an index vector that maps back to the original solution
%                   when opts.socp = 0;
%                   x = accumarray(Ech,x);
%       info.indsocp 

% How to recover the original variable x
%       after geting a solution from SeDuMi, [x;y],  for the new data Anew, bnew, cnew, Knew 
%       then, the original solution will be [accumarray(Ech,x);y]

% Author: Yang Zheng

% -------------------------------------------------------------------------
%              Input check
% -------------------------------------------------------------------------
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
        K.q = 0;
    end
    if ~isfield(opts, 'dual')
         opts.dual = 0;
    end
    if ~isfield(opts, 'keep_split')
         opts.dual = 0;
    end

% -------------------------------------------------------------------------    
%   Do not change the non-PSD part
% -------------------------------------------------------------------------
    Afree    = A(: ,1:K.f);
    Alinquad = A(:, K.f + 1:(K.l + sum(K.q)));
    %Anonpsd = A(:,1:K.f+K.l+K.q);
    %cnonpsd = c(1:K.f+K.l+K.q, :);
    cfree    = c(1:K.f, :);
    clinquad = c(K.f + 1:( K.l + sum(K.q)), :);
    
% -------------------------------------------------------------------------    
%   New SDP data
% -------------------------------------------------------------------------    
    Knew.f = K.f; Knew.l = K.l; Knew.q = K.q;
    bnew   = b;   % vector b is not channged
    Knew.s = [];
    Anew   = [];
    cnew   = [];
    Ech  = 1:K.f+K.l+K.q;      % Indexing to extract local submatrices & split cone
    Ech  = Ech(:);
    Count = K.f+K.l+K.q;  

% -------------------------------------------------------------------------    
%   The following is only used when opts.dual == 1
% -------------------------------------------------------------------------    
A_rel_free = {};

Count_rel     = 1;
Count_rel_all = 0;
i_rel_free    = [];

i_fw_free = [];
j_fw_free = [];
v_fw_free = [];

i_rel_s = [];
j_rel_s = [];
v_rel_s = [];
Count_fw_free = 0;
Count_psd = 0;

A_free_rel = cell(length(K.s), 1);
A_free_psd = cell(length(K.s), 1);


% -------------------------------------------------------------------------    
%   Approximate each PSD cone one-by-one
% -------------------------------------------------------------------------    
for PSDind = 1:length(K.s)   % multiple PSD cone

    % =====================
    % PSD data 
    % =====================
    Apsd = A(:,Count + 1:Count + K.s(PSDind)^2);  
    cpsd = c(Count + 1:Count + K.s(PSDind)^2, :);
    
    % =====================
    % get the partition 
    % =====================
    if isfield(opts, 'block')
        opts.nop = max(floor(K.s(PSDind)/opts.block), 1);
    end
    if K.s(PSDind) <= opts.nop   % the size of PSD cone must be bigger than the number of partiiton
        nop = K.s(PSDind);
    else
        nop = opts.nop;
    end

    SizeU = ceil(K.s(PSDind)/nop);
    SizeL = floor(K.s(PSDind)/nop);
    if SizeU == SizeL
        alpha = ones(nop,1)*SizeU;
    else
        x     = (K.s(PSDind) - SizeL*nop)./(SizeU-SizeL);
        alpha = [ones(x,1)*SizeU;ones(nop-x,1)*SizeL];
    end
    
    % =====================
    % split the PSD cone
    % =====================
    clique = ConeSplit(alpha);   
    
    % iterate through array
    Knew.s = [Knew.s;clique.NoElem];
    
    if opts.dual == 0 
        Ak       = cell(clique.NoC,1);
        ck       = cell(clique.NoC,1);
        for k = 1:clique.NoC   
            Tn  = cumsum([1;clique.NoElem]);
            ind = Tn(k):Tn(k+1)-1;

            Position = zeros(sum(alpha));
            Position(clique.Elem(ind),clique.Elem(ind)) = 1;
            Index  = find(Position == 1);
            Ak{k}  = Apsd(:,Index);
            ck{k}  = cpsd(Index,:);
            Anew = [Anew,Ak{k}];
            cnew = [cnew;ck{k}];
            Ech  = [Ech;Index + Count]; 
        end 
        
    elseif opts.dual == 1
        %FW*
        
        Apsd = tri_matrix(Apsd, 1);
        cpsd = tri_matrix(cpsd', 1);
        
        [i_curr, j_curr, v_curr] = find(Apsd);
        j_curr = j_curr + Count_fw_free;
        Ksi    = K.s(PSDind);
        info.ind{PSDind} = reshape(Count + tri_indexer(Ksi), [], 1);
        
        Knew.f = Knew.f + Ksi*(Ksi+1)/2;
        
        %figure out how to preallocate this
        i_fw_free = [i_fw_free; i_curr];
        j_fw_free = [j_fw_free; j_curr];
        v_fw_free = [v_fw_free; v_curr];
        
        c_fw_free_new = tri_matrix(cpsd', 1);
        
        [Mi, L] = tri_indexer(Ksi);
        Tn  = cumsum([1;clique.NoElem]);
        
        for k = 1:clique.NoC                        
            ind = Tn(k):Tn(k+1)-1;
            
            elem = clique.Elem(ind);
            
            N_curr = length(elem);
            Mi_curr = Mi(elem, elem);
            
            i_rel_free = [i_rel_free; Count_fw_free + tri_vector(Mi_curr)];
            % more work goes here
            si_curr = reshape(Count_psd + (1:N_curr^2),  N_curr, N_curr);
            %Correspondences between free variables and factor width blocks
            for i = 1:N_curr
                for j = i:N_curr
                    if i==j
                        i_rel_s = [i_rel_s; si_curr(i, j)];
                        j_rel_s = [j_rel_s; Count_rel];
                        v_rel_s = [v_rel_s; -1];                
                    else
                        i_rel_s = [i_rel_s; si_curr(i, j); si_curr(j, i )];
                        j_rel_s = [j_rel_s; Count_rel; Count_rel];
                        v_rel_s = [v_rel_s; -0.5; -0.5];                
                    end 
                    Count_rel = Count_rel + 1;
                end
            end    
            Count_psd  = Count_psd + N_curr^2;
            
        end
        Count_fw_free = Count_fw_free + Ksi*(Ksi+1)/2;
        Count_rel = Count_rel - 1;
        A_rel_free{PSDind} = sparse(i_rel_free, 1:Count_rel, ones(Count_rel, 1))';
        A_rel_psd{PSDind}  = sparse(i_rel_s, j_rel_s, v_rel_s)';
    end
    
    Count  = Count + K.s(PSDind)^2;
end 
  
Count_rel_all = Count_rel_all  + Count_rel;
  
% -------------------------------------------------------------------------    
% set up the output
% -------------------------------------------------------------------------    
  
  if opts.dual == 0
      Anew = [Afree, Alinquad, Anew];
      cnew = [cfree; clinquad; cnew];
      info.Ech = Ech;
       % Formulate an SOCP if all PSD cones are 2 by 2
        if isfield(opts,'socp') && (opts.socp == 1) && (length(find(Knew.s == 2)) == length(Knew.s))
            tic
            [Anew,bnew,cnew,Knew,indsocp] = psd2socp(Anew,bnew,cnew,Knew);
            info.socptime = toc;
            info.indsocp  = indsocp;
        end
        
  elseif opts.dual == 1
         Anew_free = sparse(i_fw_free, j_fw_free, v_fw_free);
         if opts.keep_split
            Anew = [];
            bnew = [];
            cnew = [];

            info.A_fw_free = Anew_free;
            info.A_fw_psd = sparse(length(b), sum(Knew.s.^2));
            info.A_rel_free = A_rel_free;
            info.A_rel_psd    = A_rel_psd;
            info.c_fw_free  = c_fw_free_new;
            info.c_fw_psd    = sparse(sum(Knew.s.^2), 1);
            info.new_cone = Knew;
         else
            Anew_top =  [Afree, Anew_free, Alinquad, sparse(length(b), sum(Knew.s.^2))];
            A_rel_free_diag = cell_blkdiag(A_rel_free);
            A_rel_psd_diag  = cell_blkdiag(A_rel_psd);
            A_rel = [sparse(Count_rel_all, length(cfree)) A_rel_free_diag ...
            sparse(Count_rel_all, length(clinquad)) A_rel_psd_diag];
            Anew = [Anew_top; A_rel];
            cnew =  [cfree; c_fw_free_new; clinquad; sparse( sum(Knew.s.^2), 1)];
            bnew = [b; sparse(Count_rel_all, 1)];
         end 
  end
end

