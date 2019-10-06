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
    
    Anew_lin  = [];
    Anew_quad = [];
    Anew_psd  = [];
    cnew_lin  = [];
    cnew_quad = [];
    cnew_psd  = [];
    
    Count = K.f + K.l + sum(K.q);
    
    info = cell(numPSD, 1);
    
    %set up cones array
    if ischar(cones)
        cone_str = cones;
        cones = cell(numPSD, 1);
        for i = 1:numPSD
            cones{i} = cone_str;
        end
    end
    
    for i = 1:numPSD
        K_temp = struct;
        
        
        cone_curr = cones{i};
        
        Ksi = K.s(i);
        A_curr = A(:, Count + (1:Ksi^2));
        c_curr = c(Count + (1:Ksi^2), :);
        if Ksi == 1
            Anew_lin = [Anew_lin A_curr];
            cnew_lin = [cnew_lin; c_curr];
        elseif strcmp('dd', cone_curr)
            %DD
            K_temp.s = Ksi;
            [A_dd_curr, ~, c_dd_curr, K_dd_curr, info_curr] = ...
                dd_convert(A_curr, b, c_curr, K_temp);
            
            info_curr.ind = info_curr.ind{1} + Count;
            info_curr.rays = info_curr.rays{1};
            
            Anew_lin = [Anew_lin A_dd_curr];
            cnew_lin = [cnew_lin; c_dd_curr];
            
            Knew.l = Knew.l + K_dd_curr.l;            
        elseif strcmp('sdd', cone_curr)            
            %SDD
            %SDD
            Partition = Inf;
            opts.bfw  = 1;
            opts.nop  = Partition;
            opts.socp = 1;   % second-order cones constraints

            
            K_temp.s = Ksi;
            %TODO: add indexing/recovery for sdd_info
            [A_sdd_curr, ~, c_sdd_curr, K_sdd_curr, info_curr] = ...
                factorwidth(A_curr, b, c_curr, K_temp, opts);
            
            info_curr.Ech = info_curr.Ech + Count;
            info_curr.Ech = info_curr.Ech + Count;
            
            Anew_quad = [Anew_quad A_sdd_curr];
            cnew_quad = [cnew_quad; c_sdd_curr];
            K_sdd_curr.q(K_sdd_curr.q == 0) = [];
            Knew.q = [Knew.q K_sdd_curr.q];
            
        else
            %PSD
            Anew_psd = [Anew_psd A_curr];
            cnew_psd = [cnew_psd; c_curr];
            Knew.s = [Knew.s Ksi];
            info_curr.ind = Count + (1:Ksi^2);
        end
        
        info{i} = info_curr;
        Count = Count + Ksi^2;
    end

    %output results
    Anew = [A_free_lin Anew_lin A_quad Anew_quad Anew_psd];
    cnew = [c_free_lin; cnew_lin; c_quad; cnew_quad; cnew_psd];
    bnew = b;
    
    info = [];
end

