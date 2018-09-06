function [Anew, bnew, cnew, Knew] = FactorWidth(A,b,c,K,opts)
%  Reformulating a primal semidefinte program with a block factorwidth 2
%  cone
%
%       min_{x} c^Tx
%               Ax = b
%                x \in K

% assuming K.s only, and there is only one PSD cone


% input check;
    if nargin == 5 && isfield(opts,'NoP')
        alpha = ones(opts.NoP,1)*floor(K.s/opts.NoP);
        alpha(end) = K.s - sum(alpha(1:opts.NoP-1));
    end
    if size(A,1) ~= length(b) 
        A = A';
    end


    [clique] = ConeSplit(alpha);    %% maximal cliques

    Knew.s = clique.NoElem;
    Knew.f = K.f
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
        Ak{k}  = A(:,Index);
        ck{k}  = c(Index);
        Anew = [Anew,Ak{k}];
        cnew = [cnew;ck{k}];
    end
    
    %% Formulate an SOCP if all PSD cones are 2 by 2
    if length(find(Knew.s == 2)) == length(Knew.s)
        Knew.q = ones(length(Knew.s),1)*3;
        Knew.s = 0;
        Knew.f = sum(Knew.s);
        
        for k = 1:clique.NoC   
            Tn  = cumsum([1;clique.NoElem]);
            ind = Tn(k):Tn(k+1)-1;

            Position = zeros(sum(alpha));
            Position(clique.Elem(ind),clique.Elem(ind)) = 1;
            Index  = find(Position == 1);
            Ak{k}  = A(:,Index);
            ck{k}  = c(Index);
            Anew = [Anew,Ak{k}];
            cnew = [cnew;ck{k}];
        end
    end
    
end

