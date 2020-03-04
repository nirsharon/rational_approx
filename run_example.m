% script name: "run_example"

% preparing to run
clear;
close all;

if ~exist('chebeval_scalars')
    adding_path_rational_app;
end

% function to approximate
f = @(x) abs(x-.1);
func_name = '$|x-0.1|$';

% parameters for rational app 
n = 4;
m = 4;
n_coefs = n+1;  % number of coeffs - numerator
m_coefs = m+1;  % number of coeffs - denominator
rat_name = ['rational type (',num2str(n),',',num2str(m),')'];

% parameters for optimization
LB   = 0.1;  %lower bound
UB   = 50;   % upper bound
eps1 = 1e-14;

% discretization: chebyshev points
numpt = 63;
pts   = vec(cos( pi* (2.*( numpt:-1:1) -1 ) / (2* numpt) ));

% rational app
[p,q] = RationalMinMaxOpt(f, n_coefs, m_coefs, pts, LB, UB, eps1, 0);

% prepare to plotting
LW   = 'linewidth'; 
pts  = linspace(-1,1,201).';   
func = f(pts);

% evaluate the rational approximation
p(1) = 2*p(1);
q(1) = 2*q(1);
Tp   = chebeval_scalars(p, pts ,n_coefs);
Tq   = chebeval_scalars(q, pts ,m_coefs);
rat_app = Tp(:)./Tq(:);
err_rat = func - rat_app;

% polynomial chebyshev app
coefs     = chebcoefs_app(f, n+m+1);
poly_vals = chebeval_scalars(coefs, pts ,n+m+1);
poly_name = ['Cheb. poly. deg=',num2str(n+m)];

% % polynomial interpolation
% num_int_pts = n+m;
% int_pts     = cos( pi* (2.*( num_int_pts:-1:1) -1 ) / (2* num_int_pts) );
% poly_vals   = barycentric_poly_inter2(int_pts, f(int_pts), pts);
% poly_name = ['int. poly. of deg=',num2str(n+m)];

err_poly  = func - poly_vals;
figure; 
plot(pts, err_rat,':k', LW, 3.5)
hold on
plot(pts, err_poly,'-.m', LW, 3)
ylim([-inf, max(err_poly)*1.5 ])
ylabel('Approximation error')

max_err = max(err_rat);
min_err = min(err_rat);
plot(pts(1:2:end), max_err*ones(size(pts(1:2:end))),'.r', LW, 2.8)
plot(pts, min_err*ones(size(pts)),'--r', LW, 2.8)
L = legend(rat_name, poly_name,'max error','min error','Location','southwest');
set(L,'Interpreter','latex')
set(gca,'FontSize',22)


figure; 
plot(pts, func, LW, 3)
hold on
plot(pts, rat_app,':k', LW, 3.5)
plot(pts, poly_vals,'-.m', LW, 3)
L = legend(func_name, rat_name, poly_name,'Location','NorthEast');
set(L,'Interpreter','latex')
set(gca,'FontSize',22)
