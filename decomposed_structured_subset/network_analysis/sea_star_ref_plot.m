time=[202.13; 242.40; 223.73;208.06;212.10;205.70;361.07;527.07;867.34;1683.92];

headsize = [0;1;5;10;15;20;25;30;35;40];

n_state = zeros(length(headsize), 1);
maxk_state = zeros(length(headsize), 1);

for k = 1:length(headsize)
% for k = 1:1
    curr_head = headsize(k);
    
    curr_name = strcat('sea_star_Hinf0_wide_med_', num2str(curr_head), '\\sea_star.mat');
    
    load(curr_name, 'n', 'model')
    
    n_state(k) = sum(n);
    maxk_state(k) = max(model.K.s);
end

figure
clf
plot(headsize, time/60)
xlabel('Number of agents in sea star head')
ylabel('Time (minutes)')
title('Head Size vs. H-infinity Computation Time')