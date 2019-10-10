function [Anew, bnew, cnew, Knew, info] = factorwidth(A,b,c,K,opts)
%  Reformulating a primal SDP with a block factorwidth two cone
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
%       opts.bfw     1 or 0,  block factor-width-two decomposition
%       opts.nop     integer, number of blocks in the partion alpha
%       opts.size    alternative to nop, number of entries in each block
%       opts.socp    1 or 0,  reformualte 2 by 2 PSD cone with a second-order cone
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
    Anonpsd = A(:,1:K.f+K.l+K.q);
    cnonpsd = c(1:K.f+K.l+K.q, :);
    
    %%
    Knew.f = K.f;
    Knew.l = K.l;
    Knew.q = K.q;
    bnew   = b;
    Knew.s = [];
    Anew = [];
    cnew = [];
    Ech  = 1:K.f+K.l+K.q;      % Indexing to extract local submatrices & split cone
    Ech  = Ech(:);

    
    
  Count = K.f+K.l+K.q;  
  for PSDind = 1:length(K.s)   % multiple PSD cone
      
       if isfield(opts, 'block')
            opts.nop = ceil (K.s(PSDind)/opts.block);
       end
      
       Apsd = A(:,Count + 1:Count + K.s(PSDind)^2);   % PSD data 
       cpsd = c(Count + 1:Count + K.s(PSDind)^2, :);
   
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
                x = (K.s(PSDind) - SizeL*nop)./(SizeU-SizeL);
                alpha = [ones(x,1)*SizeU;ones(nop-x,1)*SizeL];
            end
            [clique] = ConeSplit(alpha);    %% maximal cliques
            Knew.s   = [Knew.s;clique.NoElem];
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
      Count  = Count + K.s(PSDind)^2;
  end 
  
  Anew = [Anonpsd,Anew];
  cnew = [cnonpsd;cnew];
  info.Ech = Ech;
    %% Formulate an SOCP if all PSD cones are 2 by 2
    if isfield(opts,'socp') && (opts.socp == 1) && (length(find(Knew.s == 2)) == length(Knew.s))
        tic
        [Anew,bnew,cnew,Knew,indsocp] = psd2socp(Anew,bnew,cnew,Knew);
        info.socptime = toc;
        info.indsocp  = indsocp;
    end
    
end

