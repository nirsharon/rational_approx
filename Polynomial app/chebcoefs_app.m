function coefs=chebcoefs_app(func, n)
% Chebyshev fit: Given a function func, approximation such that 
%       f(x) \approx \sum_{k=0}^{n-1} c_{k}T_{k}(y) - c_{0}/2
% 
%
% Translated from Numerrical Recipes, Third edition, Section 5.8, pp. 236.
%
% Nir Sharon, October 2014.

if ~exist('n','var')
    n=50;
end
N = 2*n;
fks=zeros(N,1);
y=zeros(N,1);

coefs=zeros(n,1);


for k=0:N-1 % We evaluate the function at the n points required by (5.8.7).
    y(k+1)=cos(pi*(k+0.5)/N); 
    fks(k+1)=func(y(k+1));
end
fac=2.0/N;
for j=0:n-1 % Now evaluate (5.8.7).
    sum=0.0;
    for k=0:N-1
        sum = sum + fks(k+1)*cos(pi*j*(k+0.5)/N);
    end
    coefs(j+1)=fac*sum;
end

end