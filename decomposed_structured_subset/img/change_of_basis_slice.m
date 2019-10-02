%Plots of (decomposed) structured subsets

N = 150;
%rng(2019,'twister')
rng(200,'twister')
%rng(402, 'twister')

%Diagonal

D = [1 0 0;
     0 1 0;
     0 0 1;
     0 0 0;
     0 0 0;
     0 0 0];

%Diagonally Dominant 
 
DD_aug_C = [1  1 0  0;
      1  1 1  1;
      0  0 1  1;
      1 -1 0  0;
      0  0 1 -1; 
      0  0 0  0]/2;
  
DD_aug = [1  1 0  0 1  1;
          1  1 1  1 0  0 ;
          0  0 1  1 1  1;
          1 -1 0  0 0  0;
          0  0 1 -1 0  0;  
          0  0 0  0 1 -1]/2;  

DD = [D DD_aug];  
  
  
CDD_aug = [1   1  1  1;
           1   1  1  1;
           1   1  1  1;
           1   1 -1 -1;
           1  -1  1 -1;
           0   0  0  0]/2;
 
CDD = [DD CDD_aug];   


%PSD
R = 1;
theta = linspace(0, 2*pi, N);
phi = linspace(-pi, pi, N);

Xr = cos(phi)'*cos(theta);
X = reshape(Xr, [], 1);
Yr = cos(phi)'*sin(theta);
Y = reshape(Yr, [], 1);
Z = repmat(sin(phi)', N, 1);

PSD = [X.*X, Y.*Y,  Z.*Z,  X.*Y,  Y.*Z, X.*Z]';

%Find affine subspace through I
I = [1; 1; 1; 0; 0; 0];

A = [randn(5, 2); [0 0]];

A0 = [abs(randn(5, 1)); 0]; 
b0 = A0'*I;

%constraints and cost
K = null(A');
%c = rand(5, 1);
%c =  ones(5, 1);
%c = [1; 0; 0; 0; 0] ;
%c = [-1; 0; 1; -1; -2] ;
c = A*[ -1.5; 2.5];
%p: weights on extreme rays
c_dd = c'*DD;
At_dd = A0'*DD;
c_cdd = c'*CDD;
At_cdd = A0'*CDD;

%PSD optimization
K_psd.s = 3;
%prepare for equality constraints 
A0_psd = A0;
A0_psd(4) = A0_psd(4)/2;
A0_psd(5) = A0_psd(5)/2;
A0_psd = reshape(v2m(A0_psd), [], 1);
c_psd = reshape(v2m(c), [], 1);

[x_psd, y_psd, info_psd] = sedumi(A0_psd, b0, c_psd, K_psd);
X_psd = reshape(x_psd, 3, 3);
xv_psd = m2v(X_psd);
cost_psd = c_psd'*x_psd;

%Set up variables for change of basis
cost_dd = [];
cost_cdd = [];
DD_atom = {};
CDD_atom = {};
DD_L = {};
CDD_L1 = {};
CDD_L2 = {};
x_dd_opt = {};
x_cdd_opt = {};
ec = {};
k_max = 10;
DD_new = DD;
CDD_new = CDD;
L_prev = eye(3);
L1_prev = eye(2);
L2_prev = eye(2);

%basis changes
for k = 1:k_max
    
    %DD optimization
    n_dd = size(DD_new, 2);                                   
    K_dd.l = n_dd;

    c_dd = c'*DD_new;
    At_dd = A0'*DD_new;

    
    [p_dd, y_dd, ~] = sedumi(At_dd, b0, c_dd, K_dd);

    x_dd = DD_new*p_dd;
    X_dd = v2m(x_dd);
    cost_dd{k} = c'*x_dd;
    DD_atom{k} = DD_new;
    x_dd_opt{k} = x_dd;
    
    %new iteration of DD
    L = basis_shift(X_dd);
    L_prev = L_prev * L;
    DD_L{k} = L_prev;
    DD_new = dd_shuffle(DD, L_prev);


    %CDD optimization
    ec{k} = eig_cdd(CDD_new);
    n_cdd = size(CDD_new, 2);
    K_cdd.l = n_cdd;
    c_cdd = c'*CDD_new;
    At_cdd = A0'*CDD_new;
    [p_cdd, y_cdd, ~] = sedumi(At_cdd, b0, c_cdd, K_cdd);

    x_cdd = CDD_new*p_cdd;
    X_cdd = v2m(x_cdd);    
    
    cost_cdd{k} = c'*x_cdd;
    CDD_atom{k} = CDD_new;
    x_cdd_opt{k} = x_cdd;
    
    %new iteration of CDD
    X1 = X_cdd([1,2], [1,2]);
    X2 = X_cdd([2,3], [2,3]);

    L1 = basis_shift(X1);
    L2 = basis_shift(X2);
    
    
    %L1_prev = L1_prev * L1;
    %CDD_L1{k} = L1_prev;
    CDD_L1{k} = L1;
    %L2_prev = L2_prev * L2;
    %CDD_L2{k} = L2_prev;
    CDD_L2{k} = L2;
    CDD_new = cdd_shuffle(L1, L2);
      
end

%intersection of points in cone with affine subspace

DD_proj = {};
CDD_proj = {};
DD_opt = {};
CDD_opt = {};

kDD = {};
kCDD = {};

opt_DD = {};
opt_CDD = {};

for k = 1:k_max
    DDi = DD_atom{k};
    CDDi = CDD_atom{k};
    DD_proj0   = DDi.*b0./(A0'*DDi);
    CDD_proj0  = CDDi.*b0./(A0'*CDDi);
    DD_proj{k}  = A\DD_proj0;
    CDD_proj{k} = A\CDD_proj0;

    kDD{k} =  convhull(DD_proj{k}(1, :),    DD_proj{k}(2, :));
    kCDD{k}  = convhull(CDD_proj{k}(1, :),  CDD_proj{k}(2, :));

    
    opt_DD{k} = A\x_dd_opt{k};
    opt_CDD{k} = A\x_cdd_opt{k};
end
PSD_proj0  = PSD.*b0./(A0'*PSD);


%project onto subspace
PSD_proj = A\PSD_proj0;
I_proj = A\I;
 
opt_PSD = A\xv_psd;


%convex hulls
kPSD  = convhull(PSD_proj(1, :),  PSD_proj(2, :));

%plot
C = linspecer(5);

for k = 1:k_max
    figure(1)
    clf
    hold on

    %extreme rays
    plot(PSD_proj(1, kPSD), PSD_proj(2, kPSD),'k', 'linewidth', 2)
    plot(DD_proj{k}(1, kDD{k}), DD_proj{k}(2, kDD{k}),'color', C(1, :), 'linewidth', 2)
    plot(CDD_proj{k}(1, kCDD{k}), CDD_proj{k}(2, kCDD{k}),'color', C(2, :),'linewidth', 2)
    text(I_proj(1), I_proj(2), 'I', 'interpreter', 'latex', 'Fontsize', 20)

    %optimal points
    plot(opt_DD{k}(1), opt_DD{k}(2), '*', 'color', C(1, :), 'MarkerSize', 20, 'linewidth', 2)
    plot(opt_CDD{k}(1), opt_CDD{k}(2), '*', 'color', C(2, :), 'MarkerSize', 20,'linewidth', 2)

    plot(opt_PSD(1), opt_PSD(2), '*k', 'MarkerSize', 20, 'linewidth', 2)



    legend({'DD', 'DD(E, ?)', 'S_+'}, 'location', 'northwest', 'fontsize', 15)
    %legend({'DD', 'DD(E, ?)', 'DD_L', ' DD_L (E, ?)', 'S_+'}, 'location', 'northwest', 'fontsize', 15)
    %title('Random affine cut of structured subsets', 'Fontsize', 16)
    hold off

    xlim([min(PSD_proj(1, kPSD)),max(PSD_proj(1, kPSD))] + [-0.3 , 0.3])
    ylim([min(PSD_proj(2, kPSD)),max(PSD_proj(2, kPSD))]  + [-0.1, 0.1])

    axis square
    axis off

    keyboard;

end





