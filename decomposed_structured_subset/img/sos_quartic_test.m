sos_quartic_prep_3;

%f = x1^4 + (1+a)*x2^4 + (2+b)*x3^4 + (2-a-b)*x1^2*x2^2 / 2 + 2*b*x1*x2^3 ...
%    + (a-b)*x2^2*x3^2 + 3*a*x1*x3^3 + 2*x1^2*x2*x3;
load('quartic_process_3.mat');

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
N =  10;
%N = 4;
%N = 300;  
th = linspace(0, 2*pi, N);

%output dumps 
PSD.a = zeros(N, 1);
PSD.b = zeros(N, 1);
PSD.Q = cell(N, 1);

AGLER.a  = zeros(N, 1);
AGLER.b  = zeros(N, 1);
AGLER.Q  = cell(N, 1);

%iterate through solvers

for i = 1:N
    theta = th(i);
    model.c(1) = cos(theta);
    model.c(2) = sin(theta);

    model2.c(1) = cos(theta);
    model2.c(2) = sin(theta);

    
    %PSD
    [x,y,info] = sedumi(model.At,model.b,model.c,model.K,model.pars);
    PSD.a(i) = x(1);
    PSD.b(i) = x(2);
    K = model.K;

    PSD.Q{i} = x(K.f+1:end);

    
    %PSD(E, 0)    
    [x,y,info] = sedumi(model2.At,model2.b,model2.c,model2.K,model2.pars);
    K = model2.K;
    AGLER.a(i) = x(1);
    AGLER.b(i) = x(2);
    AGLER.Q{i} = x(K.f+1:end);
end

figure(1)
clf

%convex hulls
PSD.conv  = convhull(PSD.a, PSD.b);
AGLER.conv = convhull(AGLER.a, AGLER.b);

%plot feasibility regions
hold on 
plot(PSD.a(PSD.conv), PSD.b(PSD.conv))
plot(AGLER.a(AGLER.conv), AGLER.b(AGLER.conv))
hold off
