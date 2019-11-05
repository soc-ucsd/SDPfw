function ax = plot_region(data, color)
ax = patch('Faces',1:length(data.a),'Vertices',[data.a, data.b],...
    'EdgeColor', color, 'FaceColor', color*0.3 + [1; 1; 1]'*0.7, 'LineWidth', 2)
end