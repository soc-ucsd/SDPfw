sos_quartic_prep_3;

%f = x1^4 + (1+a)*x2^4 + (2+b)*x3^4 + (2-a-b)*x1^2*x2^2 / 2 + 2*b*x1*x2^3 ...
%    + (a-b)*x2^2*x3^2 + 3*a*x1*x3^3 + 2*x1^2*x2*x3;
%load('quartic_process_3m2.mat');
load('quartic_process_3m.mat');



%f = x1^4 + (1+a)*x2^4  + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;
%load('quartic_process0.mat');
%original
%f = x1^4 + x2^4 + a*x1^3*x2 + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;

%sparse, missing x1^3*x2 term
%f = x1^4 + (1+a)*x2^4  + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;
pars.fid = 0;
%c0 = c;

%theta = pi/6;
%number of points describing the contour
%N =  10;
%N = 30;
N = 100;
%N = 200;
%N = 300;  
th = linspace(0, 2*pi, N);


run_agler = 0;


%iterate through solvers
%output dumps


PSD = draw_feasibility(model, th);
SMALL = draw_feasibility(model_small, th);
LARGE = draw_feasibility(model_large, th);
SDD = draw_feasibility(model_sdd, th);
DD = draw_feasibility(model_dd, th);
DD_LARGE = draw_feasibility(model_large_dd, th);

if run_agler
    AGLER = draw_feasibility(model_agler, th);
    CSDD = draw_feasibility(model_csdd, th);
end

%plot feasibility regions
figure(1)
clf
C = linspecer(6);

hold on 
plot(PSD.a(PSD.conv), PSD.b(PSD.conv), 'color', 'k', 'linewidth', 2)
plot(SMALL.a(SMALL.conv), SMALL.b(SMALL.conv),  'color', C (2, :), 'linewidth', 2)
plot(LARGE.a(LARGE.conv), LARGE.b(LARGE.conv),  'color', C (3, :), 'linewidth', 2)


if  run_agler
    plot(AGLER.a(AGLER.conv), AGLER.b(AGLER.conv),  'color', C (6, :), 'linewidth', 1)
    plot(CSDD.a(CSDD.conv), CSDD.b(CSDD.conv),  'color', C (5, :), 'linewidth', 1)
    legend_str = {'S_+ (E, 0)', 'SDD(E, 0)', '[SDD, PSD]', '[PSD, SDD]', '[DD, DD]'};
else
    legend_str = {'[PSD, PSD]', '[SDD, PSD]', '[PSD, SDD]', '[SDD, SDD]', '[DD, DD]', '[PSD, DD]'};
end
plot(SDD.a(SDD.conv), SDD.b(SDD.conv),  'color', C (1, :), 'linewidth', 2)
plot(DD.a(DD.conv), DD.b(DD.conv),  'color', C (4, :), 'linewidth', 2)
plot(DD_LARGE.a(DD_LARGE.conv), DD_LARGE.b(DD_LARGE.conv),  'color', C (6, :), 'linewidth', 2)
legend(legend_str, 'location', 'northwest', 'fontsize', 15)


hold off

axis square
xlabel('a')
ylabel('b')
title('Nonnegative Feasibility Sets of f(x; a, b)', 'fontsize', 15)


function out = draw_feasibility(model, th)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.Q  = cell(N, 1);
    
    c = model.c;
    
    for i = 1:N
        theta = th(i);
        
        c(1) = cos(theta);
        c(2) = sin(theta);
        
        [x,y,info] = sedumi(model.At, model.b, c, model.K, model.pars);
        K = model.K;
        out.a(i) = x(1);
        out.b(i) = x(2);
        out.Q{i} = x(K.f+1:end);
    end
    
    out.conv = convhull(out.a, out.b);
end