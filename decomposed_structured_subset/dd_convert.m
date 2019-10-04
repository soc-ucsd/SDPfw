function [Anew, bnew, cnew, Knew, info] = dd_convert(A,b,c,K)
%  Reformulating a primal SDP with a Diagonally Dominant cone
%
%       min_{x} c^Tx
%               Ax = b
%                x \in K
%
% K can have K.f, K.l, K.q, K.s; 
%       Only replacing K.s with a block factor-width-two cone
%       and reformulating it into a standard SDP in the SeDuMi form
% Input data
%       A, b, c, K are SDP data in seudmi form
% Output data 
%       Anew, bnew, cnew, Knew, new SDP data in sedumi form
%       info.Rays   A matrix of extreme rays:
%                   
%                   x = info.Rays*x;

% How to recover the original variable x
%       after geting a solution from SeDuMi, [x;y],  for the new data Anew, bnew, cnew, Knew 
%       then, the original solution will be [info.Rays*x;y]

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
        K.q = 0;
    end

    %% Non PSD part
    A_free_lin = A(:,1:K.f+K.l);
    A_quad     = A(:, K.f+K.l + (1:K.q));
    c_free_lin = c(1:K.f+K.l);
    c_quad     = c(K.f+K.l + (1:K.q));
    
    %%setup 
    Knew.f = K.f;
    Knew.l = K.l+ sum(K.s.^2);
    Knew.q = K.q;
    bnew   = b;
    Knew.s = [];
    A_dd = [];
    c_dd = [];
    
    %%PSD part
    Count = K.f+K.l+K.q;
    Count_dd = K.f + K.l;
    info.rays = cell(length(K.s), 1);
    info.ind  = cell(length(K.s), 1);
    for PSDind = 1:length(K.s)   % multiple PSD cone
      
        Apsd = A(:,Count + 1:Count + K.s(PSDind)^2);   % PSD data 
        cpsd = c(Count + 1:Count + K.s(PSDind)^2);
         
        Ksi = K.s(PSDind);
        
        rays = dd_extreme_rays(Ksi);
        A_dd = [A_dd Apsd*rays];
        c_dd = [c_dd (cpsd'*rays)'];
        
        info.ind{PSDind} = Count_dd + (1:Ksi^2);
        info.rays{PSDind} = rays;
        Count_dd = Count_dd + Ksi^2;
        Count = Count + Ksi^2;
    end
    
    Anew = [A_free_lin A_dd A_quad];
    cnew = [c_free_lin; c_dd; c_quad];

end