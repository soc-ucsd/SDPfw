c = linspecer(2);

figure(340)
clf
hold on
rx = [0,1,1,0,0];
ry = [0, 0,-1,-1,0];
rect = [rx; ry];

%coarse decomposition
coarse = [30, 20, 40, 20, 40];
accum_len = 0;
for i = 1:length(coarse)
   curr_len = coarse(i);
   curr_rect = curr_len*rect + [1;-1]*accum_len;
   patch(curr_rect(1, :), curr_rect(2, :), 'm', 'EdgeColor', 'none')
   accum_len = accum_len+curr_len;
end

%fine decomposition
for i = 1:15
    curr_rect = 10*rect + (i-1)*[1;-1]*10;
   patch(curr_rect(1, :), curr_rect(2, :), 'b', 'EdgeColor', 'none')
end
rect_arrow_bot = [160*rx; -150 + 10*ry];
patch(rect_arrow_bot(1, :), rect_arrow_bot(2, :), 'b', 'EdgeColor', 'none')

rect_arrow_side= [150 + 10*rx; 160*ry];
patch(rect_arrow_side(1, :), rect_arrow_side(2, :), 'b', 'EdgeColor', 'none')

axis square
axis off