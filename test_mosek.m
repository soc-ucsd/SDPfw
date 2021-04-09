%%  test the usage of factorwitdth.m using Mosek

clear;close all
Partition = [2:2:10];         % partition of blocks

file = 'SedumiDataEx10';
load(['SeDuMiData\', file '.mat']);

Tcon = zeros(length(Partition),1);
Time = zeros(length(Partition),1);
Cost = zeros(length(Partition),1);

for k = 1:length(Partition)
    opts.bfw = 1;
    opts.nop = Partition(k);
    opts.socp = 1;
    
    % reformulating the SDP
    Ts = tic;
    [Anew, bnew, cnew, Knew, infofw] = factorwidth(A',b,c,K,opts);
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
subplot(1,2,1)
plot(Partition, Time,'*','markersize',10);hold on; plot(Partition, Time,'b','linewidth',1.2);
xlabel('Number of partition','interpreter','latex')
ylabel('Solver time consumption (s)','interpreter','latex')

subplot(1,2,2)
plot(Partition, Cost,'*','markersize',10);hold on; plot(Partition, Cost,'b','linewidth',1.2);
xlabel('Number of partition','interpreter','latex')
ylabel('Cost value','interpreter','latex')


