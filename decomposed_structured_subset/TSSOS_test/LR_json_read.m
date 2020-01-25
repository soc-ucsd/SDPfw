data = jsondecode(fileread('LR12_trans_box.json'));

N = data.n;

fbasis = reshape(data.fbasis, N, []);
gbasis = cell(N, 1);
for i = 1:N
    gb_curr = data.gbasis(i, :);
    gbasis{i} = reshape(gb_curr, N, []);
end


fblocks = data.fblocks;
gblocks = data.gblocks;

%formulas 22 and 23 in the TSSOS paper in order to get the Psatz in YALMIP
%each block in fblocks indexes into the support set in fbasis
%likewise for gblocks
