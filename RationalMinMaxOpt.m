function [p, q, z] = RationalMinMaxOpt(f, n, m, pts, LB, UB, prc, vrb)
% Calculating the uniform best rational approx of type (n,m)
% via optimization with deviation precision 'eps'
%
% Input:
%   f   - the function to be approximated. A function handler
%   n,m - the rational approx parameters = maximum degree (numer.,deno.)
%   pts - discretization points
%   LB  - lower bound on the denominator (away from zero)
%   UB  - upper bound on the denominator
% 	prc - precision of the bisection (maximum deviation accuracy)
%   vrb - flag for verbose run
% Ouput:
%   p,q - the rational approx coefficients 
%   z   - the maximal deviation
% Run example:
%  [p, q, z] = RationalMinMaxOpt(@(x) abs(x), 4, 4, linspace(-1,1,20), .1, 100, 1e-10, 1)
%
% Elior Kalfon, Nir Sharon, Feb 2020

% verbose run
if nargin < 8
    vrb = 0;
end
% deviation precision
if nargin < 7
    prc = 10^-14;
end
% upper bound on denominator (not necessary in general)
if nargin < 6
    UB = 1000*LB;
end

% "Chebyshev" Vandermonde matrix
I = eye(n); I(1) = 2;
Tn = zeros(numel(pts),n);
for deg=0:(n-1)
    Tn(:,deg+1) = chebeval_scalars(I(deg+1,:), pts, deg+1);
end
I = eye(m); I(1) = 2;
Tm = zeros(numel(pts),m);
for deg=0:(m-1)
    Tm(:,deg+1) = chebeval_scalars(I(deg+1,:), pts, deg+1);
end

% lower bound
uL  = 0;

% upper bound by polynomial interpolation (degree n, chebyshev pts)
int_pts = vec(cos( pi* (2.*( n:-1:1) -1 ) / (2*n) ));
uH      = max(abs( f(pts) - barycentric_poly_inter2(int_pts, f(int_pts), pts)));

while(uH-uL)>prc
    z=(uH+uL)/2;
    %If the optimal result is lower than distance z 
    if(checkVal(f,z,n,m,pts,Tn,Tm,LB,UB,vrb))
        uH=z;
    else
        uL=z;
    end
end
%Calculates the optimal p,q 
[p,q,~]=LpRat(f,uH,n,m,pts,Tn,Tm,LB,UB,vrb);
end

function bool=checkVal(f,z,n,m,T,Tn,Tm,DLB,DUB,vrb)
    [~,~,u,exitflag]=LpRat(f,z,n,m,T,Tn,Tm,DLB,DUB,vrb);
    bool= (u<= 1e-12);% (u<= 10^-15); % exitflag; % To be refined.
end

function [p,q,u,exitflag] = LpRat(f,z,n,m,T,Tn,Tm,DLB,DUB,vrb)
%cond1 and cond2 are the one by column multiplication of each site and a row of
%chebyshev matrix
cond1=multRow(f(T)+z,Tm);
cond2=multRow(f(T)-z,Tm);
%A is a matrix of constraints corresponds to the inequation Ax<b
A=[Tn -cond1 -ones(length(T),1);
    -Tn cond2 -ones(length(T),1);
    zeros(size(Tn)) -Tm zeros(length(T),1);
    zeros(size(Tn)) Tm zeros(length(T),1)];

b=[zeros(2*size(Tn,1),1);-DLB*ones(length(T),1);
    DUB*ones(length(T),1)];

% lb=[]; %[-inf.*ones(1,n+m+1)];
lb=[-inf.*ones(1,n) 0 -inf.*ones(1,m-1) -1];

%ub=[inf.*ones(1,n) 1 inf.*ones(1,m)];
ub=[inf.*ones(1,n) inf.*ones(1,m+1)];

%A dummy varibale in order to arrange in one objective function
%obj=[zeros(1,n+m) 1];
obj=[zeros(1,n+m) 1];

if vrb
    options = optimoptions('linprog','Display','iter');
else
    options = optimoptions('linprog','Display','none');
end
if vrb
    [x,~,exitflag,output] = linprog(obj,A,b,[],[],lb,ub,options);
    output
else
    [x,~,exitflag,~] = linprog(obj,A,b,[],[],lb,ub,options);
end
p=x(1:n);
q=x(n+1:m+n);
u=x(n+m+1);

end

function M=multRow(F,Tm)
M=zeros(size(Tm));
for i=1:length(F)
    M(i,:)=F(i).*Tm(i,:);
end
end
