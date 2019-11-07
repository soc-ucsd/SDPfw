% head = 60;      %size of central 'head'
% knuckle = 10;    %size of each knuckle
% t = 4;          %#links between head and first knuckle
% t_k = 2;        %#links between subsequent knuckles
% N = 6;          %#arms
% k = 6;          %#knuckles per arm

%load('sea_star_hinf.mat')
%load('sea_star_Hinf.mat')
%load('sea_star_Hinf_small0.mat')
%load('sea_star_Hinf_medium1.mat')

load('sea_star_H2_small1.mat')

%load('sea_star_Hinf_medium.mat')
Js = LOP.J.s';
cone = cell(length(Js), 1);

% %Cone thresholds
%[0, 22, 70]
% cone_thresh = 70;

cone_psd = cone_list(Js, 0, 'psd');
cone_dd  = cone_list(Js, 0, 'dd');
cone_22_dd  = cone_list(Js, 22, 'dd');
cone_50_dd  = cone_list(Js, 50, 'dd');
%cone_70_dd  = cone_list(Js, 70, 'dd');
cone_100_dd  = cone_list(Js, 100, 'dd');

use_mosek = 1;

[Hout, res, time_solve, time_convert] = run_model_star(LOP, cone_psd, use_mosek);
[Hout_dd, res_dd, time_solve_dd, time_convert_dd] = run_model_star(LOP, cone_dd, use_mosek);
[Hout_dd_22, res_dd_22, time_solve_dd_22, time_convert_dd_22] = run_model_star(LOP, cone_22_dd, use_mosek);
[Hout_dd_50, res_dd_50, time_solve_dd_50, time_convert_dd_50] = run_model_star(LOP, cone_50_dd, use_mosek);
%[Hout_dd_70, res_dd_70, time_solve_dd_70, time_convert_dd_70] = run_model_star(LOP, cone_70_dd, use_mosek);
[Hout_dd_100, res_dd_100, time_solve_dd_100, time_convert_dd_100] = run_model_star(LOP, cone_100_dd, use_mosek);

output = [Hout_dd Hout_dd_22 Hout_dd_50 Hout_dd_100 Hout]
%output = [Hout_dd Hout_dd_22 Hout_dd_50 Hout_dd_70 Hout_dd_100 Hout]
%output = [Hout_dd Hout_dd_22 Hout_dd_100]





