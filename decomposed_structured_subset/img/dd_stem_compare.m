x1 = 3;
x3 = 1;
x4 = 4;


x2_dd  = x1 + x3 + x4;
x2_cdd = max(x1, x3+x4);

sz = 5000;

sz_stem = 50;
fz = 26;
r = 1:4;

figure(2)
clf
hold on

stem([1,3,4], [x1, x3, x4], '.k', 'filled', 'markersize', sz_stem) 

d = 0.1;
text(1.1, x1, strcat('$x_{12}=', num2str(x1), '$'), 'Interpreter', 'Latex', 'Fontsize', fz)
text(3.1, x3, strcat('$x_{23}=', num2str(x3), '$'), 'Interpreter', 'Latex', 'Fontsize', fz)
text(3.45, x4, strcat('$x_{24}=', num2str(x4), '$'), 'Interpreter', 'Latex', 'Fontsize', fz)

ymax = x2_dd + 4;
lw = 3;

line([2, 2], [x2_cdd, ymax], 'color',C(2, :), 'linewidth', lw)
scatter(2, x2_cdd, sz, C(2, :), '.')
text(2.1, x2_cdd, strcat('$x_{22}\geq', num2str(x2_cdd), '$'), 'Interpreter', 'Latex', 'Fontsize', fz,  'color',C(2, :))
text(2.1, x2_cdd+1, strcat('$X \in \{\mathcal{D}\mathcal{D}^2, \mathcal{D}\mathcal{D}^3\}$'), 'Interpreter', 'Latex', 'Fontsize', fz,  'color',C(2, :))

line([2, 2], [x2_dd, ymax], 'color', C(1, :), 'linewidth', lw)
scatter(2, x2_dd, sz, C(1, :), '.')
text(2.1, x2_dd, strcat('$x_{22}\geq', num2str(x2_dd), '$'), 'Interpreter', 'Latex', 'Fontsize', fz,  'color',C(1, :))
text(2.1, x2_dd+1, strcat('$X \in \mathcal{D}\mathcal{D}^4$'), 'Interpreter', 'Latex', 'Fontsize', fz,  'color',C(1, :))
hold off

axis off