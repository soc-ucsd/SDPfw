function [Anew, bnew, cnew, Knew] = FactorWidth(A,b,c,K,opts)
%  Reformulating a primal semidefinte program with a block factorwidth 2
%  cone
%
%       min_{x} c^Tx
%               Ax = b
%                x \in K

% assuming K.s only, and there is only one PSD cone
% can have K.f, K.l, K.q


% input check;
    if nargin == 5 && isfield(opts,'NoP')
        alpha = ones(opts.NoP,1)*floor(K.s/opts.NoP);
        alpha(end) = K.s - sum(alpha(1:opts.NoP-1));
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
    Anonpsd = A(:,1:K.f+K.l+K.q);
    Apsd    = A(:,K.f+K.l+K.q+1:end);
    
    cnonpsd = c(1:K.f+K.l+K.q);
    cpsd    = c(K.f+K.l+K.q+1:end);
    %%

    [clique] = ConeSplit(alpha);    %% maximal cliques

    Knew.s = clique.NoElem;
    Knew.f = K.f;
    Knew.l = K.l;
    Knew.q = K.q;
    bnew = b;
    Ak   = cell(clique.NoC,1);
    ck   = cell(clique.NoC,1);
    Anew = [];
    cnew = [];
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
    end
    
    Anew = [Anonpsd,Anew];
    cnew = [cnonpsd;cnew];
    
    
    %% Formulate an SOCP if all PSD cones are 2 by 2
    if isfield(opts,'socp')&&(opts.socp == 1) && (length(find(Knew.s == 2)) == length(Knew.s))
        [data,Knew] = sedumi2scs(Anew',bnew,cnew,Knew);
        
        Asocp = data.A';
        bsocp = -data.c;
        csocp = data.b;
        
        H = [1 0 1;0 2 0; 1 0 -1]^(-1);
        for k = 1:clique.NoC
            ind           = 3*(k-1)+1:3*k;    % the socp cone size is 3
            Asocp(:,ind)  = Asocp(:,ind)*H;
            csocp(ind)    = H'*csocp(ind);
        end
               
        Knew.q = ones(1,length(Knew.s))*3;
        Knew.s = 0;
        Anew = Asocp;
        bnew = bsocp;
        cnew = csocp;
        
        %% add nonnegative constraints
    end
    
end

