function [Anew, bnew, cnew, Knew, info] = dd_star_convert(A,b,c,K, keep_split)
%  Reformulating a primal SDP with a Diagonally Dominant * cone
%  The DD* cone is larger than the PSD cone, and will have a lower
%  objective value (is an outer approximation).
%
%       min_{x} c^Tx
%               Ax = b
%                x \in K*
%
% K can have K.f, K.l, K.q, K.s; 
% Input data
%       A, b, c, K are SDP data in seudmi form
% Output data 
%       Anew, bnew, cnew, Knew, new SDP data in sedumi form
%       info for recovery 
%
% keep_split = 0: perform decomposition and new SDP
% keep_split = 1: don't merge together relation SDPs

%% Input check

    if nargin <= 4
        keep_split = 0;
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
        K.q = 0;
    end

    %% Non PSD part
    A_free = A(:,1:K.f);
    A_lin  =  A(:, K.f + (1:K.l));
    A_quad = A(:, K.f+K.l + (1:K.q));
    c_free = c(1:K.f);
    c_lin  = c(K.f + (1:K.l));
    c_quad = c(K.f+K.l + (1:K.q));
    
    %%setup 
    Knew.f = K.f + sum(K.s.*(K.s+1)/2);
    Knew.l = K.l+ sum(K.s.^2);
    Knew.q = K.q;
    bnew   = b;
    Knew.s = [];
    %A_dd = [];
    %c_dd = [];
    %A_dd_free = [];
    i_dd_free = [];
    j_dd_free = [];
    v_dd_free = [];
    c_dd_free = [];
    A_dd_lin  = sparse(length(b), sum(K.s.^2));
    c_dd_lin  = sparse(sum(K.s.^2), 1);
    
    A_rel_free = cell(length(K.s), 1);
    A_rel_lin  = cell(length(K.s), 1);
    
    %%PSD part
    Count = K.f+K.l+sum(K.q);
    %Count_dd
    Count_dd_free = 0;
%    Count_dd_lin  = K.f + K.l;
    %Count_dd_lin  = 0;
    
    %info.rays = cell(length(K.s), 1);
    info.ind  = cell(length(K.s), 1);
    for PSDind = 1:length(K.s)   % multiple PSD cone
        Ksi = K.s(PSDind);
        num_dd = Ksi^2;
        
        Apsd = A(:,Count + 1:Count + Ksi^2);   % PSD data
        cpsd = c(Count + 1:Count + Ksi^2, :);
        
        %[At_psd_dd, c_dd_free_new] = svecData_mod(Apsd', cpsd, Ksi);
        %rays_svec, ~] = svecData_mod(rays', cpsd, Ksi);
        At_psd_dd = tri_matrix(Apsd, 1);
        c_dd_free_new = tri_matrix(cpsd', 1);
        
        %extreme rays for DD cone
        rays = dd_extreme_rays(Ksi);
        %rays(:, Ksi+3:end) = 0;
        %raysT  = tri_matrix(rays', 1);
        raysT  = tri_matrix(rays', 1);
        %raysT = zeros(size(raysT));
        
        num_svec = length(c_dd_free);
        
        %cost, and standard affine constraints <Ai, X> = bi
        %A_dd_free = [A_dd_free At_psd_dd];
        [i_curr, j_curr, v_curr] = find(At_psd_dd);

        j_curr = j_curr + Count_dd_free;
        
        info.ind{PSDind} = reshape(Count + tri_indexer(Ksi), [], 1);
        
        Count_dd_free = Count_dd_free + Ksi*(Ksi+1)/2;
        %figure out how to preallocate this
        i_dd_free = [i_dd_free; i_curr];
        j_dd_free = [j_dd_free; j_curr];
        v_dd_free = [v_dd_free; v_curr];
        
        
        c_dd_free = [c_dd_free; c_dd_free_new'];
        
        %A_dd_lin  = [A_dd_lin; rays_svec'];
        %c_dd_lin  = [c_dd_lin; sparse(num_dd, 1)];                                
        %A_dd_lin  = [A_dd_lin sparse(length(b), num_dd)];
        %now relate the new variables together, rays*free = lin
        
        A_rel_free{PSDind} = raysT;
        A_rel_lin{PSDind} = -speye(num_dd);
        
        %A_dd_free = 0;
        
        %A_dd = [A_dd Apsd*rays];
        %c_dd = [c_dd; rays*cpsd];
        
        %info.ind{PSDind} = Count_dd + (1:Ksi^2);
        %info.rays{PSDind} = rays;
        
        %Count_dd_free = Count_dd_free + num_svec;
        %Count_dd_lin  = Count_dd_lin + num_svec + num_dd;
        
        %Count_dd = Count_dd + Ksi^2;
        Count = Count + num_dd;        
    end
    
    A_dd_free = sparse(i_dd_free, j_dd_free, v_dd_free, length(b), sum(K.s .*(K.s + 1)/2));
    
    if keep_split
        Anew = [];
        bnew = [];
        cnew = [];
        info.A_dd_free = A_dd_free;
        info.A_dd_lin  = A_dd_lin;
        info.A_rel_free  = A_rel_free;
        info.A_rel_lin = A_rel_lin;
        info.c_dd_free = c_dd_free;
        info.c_dd_lin  = c_dd_lin;
    else        
        Anew = [A_free A_dd_free A_lin  A_dd_lin A_quad];
        cnew = [c_free; c_dd_free; c_lin; c_dd_lin; c_quad];


        A_rel_free_diag = cell_blkdiag(A_rel_free);
        A_rel_lin_diag  = cell_blkdiag(A_rel_lin);

        num_rel = size(A_rel_free_diag, 1);
        A_rel = [sparse(num_rel, length(c_free)) A_rel_free_diag ...
                 sparse(num_rel, length(c_lin))  A_rel_lin_diag ...
                 sparse(num_rel, length(c_quad))];
        b_rel = sparse(num_rel, 1);

        Anew = [Anew; A_rel];
        bnew = [bnew; b_rel];
    end
    
    %now for the relations between entries
    
end