%load('block_arrow_change.mat')

figure(1)
clf
C = linspecer(5);
hold on
plot(1:iter_max, Cost_Matrix(:, 1)/1e4, '-s','color', C(1, :), 'linewidth', 2, 'markersize', 10)
plot(1:iter_max, Cost_Matrix(:, 2)/1e4, '-o', 'color', C(2, :), 'linewidth', 2, 'markersize', 10)
plot(1:iter_max, Cost_Matrix(:, 3)/1e4, '-*',  'color', C(3, :), 'linewidth', 2, 'markersize', 10)

plot(xlim, [cost0, cost0]/1e4, 'k--', 'linewidth', 2)
hold off
legend_list = {'$B_5$', '$B_5(E_F, ?)$', '$B_5(E, ?)$', '$S_+$'};
legend(legend_list, 'Interpreter', 'latex', 'fontsize', 18)
title('Decomposed Change of Basis Comparision', 'fontsize', 18)

xlabel('Iteration')
ylabel('Cost')