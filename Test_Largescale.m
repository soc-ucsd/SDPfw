%%
clc;clear
% load TestEx

fileName{1} = 'rs365';
fileName{2} = 'rs1555';
fileName{3} = 'rs1907';
fileName{4} = 'maxG11';
fileName{5} = 'maxG32';

for Index = 4:5
    file = fileName{Index};
    load([file '.mat']);

Partition = 10:10:50;
Partition = [2,4,6,8,Partition];
Tcon = zeros(length(Partition),1);
Time = zeros(length(Partition),1);
Cost = zeros(length(Partition),1);

for k = 1:length(Partition)
    opts.NoP = Partition(k);
    Ts = tic;
    [Anew, bnew, cnew, Knew] = FactorWidth(At,b,c,K,opts);
    Tcon(k) = toc(Ts);

    prob1          = convert_sedumi2mosek(Anew', bnew, cnew,Knew); 
    try
        Tmosek         = tic;
        [rcode1, res1] = mosekopt('minimize info', prob1);
        Tmosek1        = toc(Tmosek);
        Time(k)  = Tmosek1;
        Cost(k)  = res1.sol.itr.pobjval;
    catch
        
    end
    
    save(['Result_' file],'Time','Cost','Tcon')
end

    
end