% =========================================================================
%      test factorwitdth.m using sedumi
% =========================================================================


clear; close all;

% load SDP data
file = 'SedumiDataEx10';
load(['SeDuMiData\', file '.mat']);

% Set partitions
Partition = [2:2:10];                  % partition of blocks
Tcon = zeros(length(Partition),1);    % time for conversion
Time = zeros(length(Partition),1);    % time for sedumi
Cost = zeros(length(Partition),1);    % cost value

x = zeros(length(c),length(Partition));
y = zeros(length(b),length(Partition));  % record original variables

% run the tests
for k = 1:length(Partition)
    opts.bfw  = 1;
    opts.nop  = Partition(k);
    opts.socp = 1;   % second-order cone constraints
    
    % Approximate the SDP using a block factor-width cone
    % Reformualte the block factor-width cone program into a standard SDP
    Ts      = tic;
    [Anew, bnew, cnew, Knew, infofw] = factorwidth(A,b,c,K,opts);
    Tcon(k) = toc(Ts);
    
    % solve the new SDP using sedumi
    Tsedumi      = tic;
    [xn,yn,info] = sedumi(Anew,bnew,cnew,Knew);
    Time(k)      = info.cpusec;
    Cost(k)      = cnew'*xn;

    % variables corresponds to the original SDP data, A, b, c, K
    % y(:,k) = yn;
    % x(:,k) = accumarray(infofw.Ech,xn);
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

