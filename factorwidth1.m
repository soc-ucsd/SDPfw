function [Anew, bnew, cnew, Knew, Ech] = factorwidth1(A,b,c,K,opts)
%  Reformulating a primal SDP with a block factorwidth two cone
%
%       min_{x} c^Tx
%               Ax = b
%                x \in K
%
% K can have K.f, K.l, K.q, K.s
% Only replacing K.s with a block factor-width-two cone
% Reformulating it into a standard SDP in the SeDuMi form

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
    cnonpsd = c(1:K.f+K.l+K.q);
    
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
      
       Apsd = A(:,Count + 1:Count + K.s(PSDind)^2);   % PSD data 
       cpsd = c(Count + 1:Count + K.s(PSDind)^2);
   
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
                ck{k}  = cpsd(Index);
                Anew = [Anew,Ak{k}];
                cnew = [cnew;ck{k}];
                Ech  = [Ech;Index + Count];
            end
      Count  = Count + K.s(PSDind)^2;
  end 
  
  Anew = [Anonpsd,Anew];
  cnew = [cnonpsd;cnew];
    
    %% Formulate an SOCP if all PSD cones are 2 by 2
    if isfield(opts,'socp') && (opts.socp == 1) && (length(find(Knew.s == 2)) == length(Knew.s))
        % to do
    end
    
end

