%%  test the usage of factorwitdth.m using Mosek

clc;clear
Partition = [2,5,10,20,50,100];         % partition of blocks

file = 'SedumiDataEx15';
load(['SeDuMiData\', file '.mat']);

Tcon = zeros(length(Partition),1);
Time = zeros(length(Partition),1);
Cost = zeros(length(Partition),1);

for k = 1:length(Partition)
    opts.bfw = 1;
    opts.nop = Partition(k);
    
    % reformulating the SDP
    Ts = tic;
    [Anew, bnew, cnew, Knew, Ech] = factorwidth(A',b,c,K,opts);
    Tcon(k) = toc(Ts);
    
    % solve the new SDP using Mosek
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
end

figure
plot(Partition, Time,'*','markersize',10);hold on; plot(Partition, Time,'b','linewidth',1.2);
xlabel('Number of partition','interpreter','latex')
ylabel('Solver time consumption (s)','interpreter','latex')

