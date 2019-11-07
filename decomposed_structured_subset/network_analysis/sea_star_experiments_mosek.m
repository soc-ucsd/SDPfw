%load('sea_star_H2_small.mat')
%load('sea_star_Hinf0_small.mat')
%load('sea_star_Hinf0_medium.mat')

%fname = 'sea_star_H2_tiny.mat';
fname = 'sea_star_Hinf0_tiny.mat';
%fname = 'sea_star_H2_small.mat';
%fname = 'sea_star_Hinf0_medium.mat';

[filepath,name,ext] = fileparts(fname);

load(fname);

Js = LOP.J.s';


cones = {'dd', 1, 2, 4, 7};
thresh = [0, 4, 10, 17];

%cones = {'dd', 'sdd', 3, 6, 12, 18, Inf};
%cones = {'dd', 6, 18, Inf};
%thresh = [0, 22, 50, Inf];
%thresh = [0, 11, 35, Inf];
%thresh = [0, 11, 30, 50,  Inf];
%thresh = [0, 11, 33,  56,  Inf];
%thresh = [0, 11, 33,  56];

Ncones = length(cones);
Nthresh = length(thresh);

CONE = cell(Ncones, Nthresh);
RES  = cell(Ncones, Nthresh);
for i = 1:Ncones
    for j = 1:Nthresh
        CONE{i,j} = struct;
        CONE{i,j}.cone = cone_list(Js, thresh(j), cones{i});
    end
end

cone = cell(length(Js), 1);

% %Cone thresholds
%[0, 22, 70]
% cone_thresh = 70;

% cone_psd = cone_list(Js, 0, 'psd');
% cone_dd  = cone_list(Js, 0, 'dd');
% cone_22_dd  = cone_list(Js, 22, 'dd');
% cone_50_dd  = cone_list(Js, 50, 'dd');
% %cone_70_dd  = cone_list(Js, 70, 'dd');
% cone_100_dd  = cone_list(Js, 100, 'dd');

use_mosek = 1;

output = NaN*ones(Ncones, Nthresh);


[CONE0.Hout, RES0, CONE0.time_solve, CONE0.time_convert]...
            = run_model_star(LOP, 'psd', use_mosek);
fprintf('Cone: PSD \t Hinf: %3f\n', CONE0.Hout)

for i = 1:Ncones
    for j = 1:Nthresh
        [CONE{i,j}.Hout, RES{i,j}, CONE{i,j}.time_solve, CONE{i,j}.time_convert]...
            = run_model_star(LOP, CONE{i,j}.cone, use_mosek);              
             
        output(i, j) = CONE{i,j}.Hout;
              %output
             %{cones{i}, thresh(j), output(i,j)}
         fprintf('Cone: %s \t Thresh:  %d \t Hinf: %3f\n', num2str(cones{i}), thresh (j), output(i,j))
%         else
%             CONE{i,j}.Hout = NaN;        
%         end
    end
end

[output ones(Ncones, 1)*CONE0.Hout]

outname = strcat(filepath,'output_',name,ext);
save(outname, 'CONE', 'CONE0', 'cones', 'thresh')





