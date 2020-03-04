function fv=chebeval_scalars(coef, pts ,m)
% Chebyshev evaluation: The Chebyshev polynomial
%       \sum_{k=0}^{m-1} c_{k}T_{k}(x) - c_{0}/2 , x\in pts
%
% Translated from Numerical Recipes, Third edition, Section 5.8, pp. 237.
%
% NS, Dec 19. Based on a joint work with Y.Shkolnisky


if max(abs(pts))>1
    error('Numbers outside the segment [-1,1]');
end

if m>numel(coef) % Lowest order of apprixmation is one, which corrponds to the constant approximation
    error('Approximation order is too high for the precomputed C');
end

if m<1
    error('Approximation order must be greater than 1');
end

d  = zeros(size(pts));
dd = zeros(size(pts));
n  = numel(pts,1);
y2 = 2*pts; 

for j=m-1:-1:1 % Clenshaw's recurrence.
    sv = d;
    d  = y2.*d-dd+coef(j+1);
    dd = sv;
end
fv= pts.*d - dd + 0.5*coef(1);

end
