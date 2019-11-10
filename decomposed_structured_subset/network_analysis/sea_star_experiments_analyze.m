%load('output_sea_star_Hinf0_large.mat')
load('output_sea_star_Hinf0_verylarge.mat')

Ncones = length(cones);
Nthresh = length(thresh);

Hout = zeros(Ncones, Nthresh);
time_solve = zeros(Ncones, Nthresh);
time_convert = zeros(Ncones, Nthresh);


for i = 1:Ncones
    for j = 1:Nthresh
        Hout(i, j) = CONE{i,j}.Hout;
        time_solve(i, j) = CONE{i,j}.time_solve;
        time_convert(i, j) = CONE{i,j}.time_convert;               
    end
end

Hout0 = CONE0.Hout;
time_solve0 = CONE0.time_solve;
time_convert0 = CONE0.time_convert;

time = time_solve + time_convert;
time0 = time_solve0 + time_convert0;

%cost_table= latex(vpa(sym(time_solve),4));

total_time = sum(sum(time)) + time0;

cost_table = latex(vpa(sym(Hout),4));
time_table = latex(vpa(sym(time),4));

min(min(time))/time0;