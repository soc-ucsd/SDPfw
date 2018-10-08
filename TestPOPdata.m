%%
clc;clear

folder = 'SeDuMiData\';

fileName{1} = 'SedumiDataEx10';
fileName{2} = 'SedumiDataEx15';
fileName{3} = 'SedumiDataEx20';
fileName{4} = 'SedumiDataEx25';
fileName{5} = 'SedumiDataEx30';
fileName{6} = 'SedumiDataEx35';
fileName{7} = 'SedumiDataEx40';

Partition = [2,4,6,8,10,20,30,40,50,75,100];

for Index = 1:length(fileName)
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
