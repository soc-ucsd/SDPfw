fileName{1} = 'NewResult_SedumiDataEx10';
fileName{2} = 'NewResult_SedumiDataEx15';
fileName{3} = 'NewResult_SedumiDataEx20';
fileName{4} = 'NewResult_SedumiDataEx25';
fileName{5} = 'NewResult_SedumiDataEx30';
% fileName{6} = 'NewResult_SedumiDataEx35';
% fileName{7} = 'NewResult_SedumiDataEx40';

Partition = [2,4,6,8,10,20,30,40,50,75,100];
folder = 'SeDuMiData\';

Time_all = zeros(length(Partition)-1, 5);
Cost_all = zeros(length(Partition)-1, 5);
figure(1)
clf 
hold on
figure(2)
clf 
hold on
inds_= [1:3, 5:11];
for ii = 1: 5
    file = fileName{ii};
    load([folder, file '.mat']);
    Time_all(:,ii) = Time(inds_);
    Cost_all(:,ii) = Cost(inds_);
    figure(1); plot(Partition(inds_), log10(Time(inds_)))
    figure(2); plot(Partition(inds_), (Cost(inds_)))
end