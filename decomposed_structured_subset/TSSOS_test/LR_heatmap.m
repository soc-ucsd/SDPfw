%SC = load('LR120_output_box_1_2.mat', 'CONE');
%CONE_cons = SC.CONE;

cost_cons = cellfun(@(x) x.cost, CONE_cons);
time_cons = cellfun(@(x) x.time_solve + x.time_convert, CONE_cons);

%UC = load('LR120_output_uncons.mat', 'CONE');
%CONE_uncons = UC.CONE;

cost_uncons = cellfun(@(x) x.cost, CONE_cons);
time_uncons = cellfun(@(x) x.time_solve + x.time_convert, CONE_cons);