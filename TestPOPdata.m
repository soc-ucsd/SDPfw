%%
clc;clear

folder = 'SeDuMiData\';

fileName{1} = 'sedumiPOPbt_10';
fileName{2} = 'sedumiPOPbt_15';
fileName{3} = 'sedumiPOPbt_20';
fileName{4} = 'sedumiPOPbt_30';
fileName{5} = 'sedumiPOPbt_40';

Partition = [2,4,6,8,10,20,30,40,50];

for Index = 1:5
    file = fileName{Index};
    load([folder, file '.mat']);

    Tcon = zeros(length(Partition),1);
    Time = zeros(length(Partition),1);
    Cost = zeros(length(Partition),1);

    for k = 1:length(Partition)
        opts.NoP = Partition(k);
        Ts = tic;
        [Anew, bnew, cnew, Knew] = FactorWidth(A',b,c,K,opts);
        Tcon(k) = toc(Ts);

        prob1          = convert_sedumi2mosek(Anew', bnew, cnew,Knew); 
        try
            Tmosek         = tic;
            [rcode1, res1] = mosekopt('minimize info', prob1);
            Tmosek1        = toc(Tmosek);
            Time(k)  = Tmosek1;
            Cost(k)  = res1.sol.itr.pobjval;
        catch
            Time(k)  = NaN;
            Cost(k)  = NaN;
        end

        save([folder, 'Result_' file],'Time','Cost','Tcon')
    end
end
