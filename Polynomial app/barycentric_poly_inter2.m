function [vals] = barycentric_poly_inter2(xj, yj, x)
% Polynomial interpolation using barycentric representation
%
% Input:
%   x   -- the interpolation points
%   y -- the values of the function over interpolation pts
%   pts  -- the new points to evaluate the interpolant
% Output:
%   vals    -- the interpolation polynomial values over "eval_pts"
% NS, Jan 19

n   = length(xj);
thd = 1e-10;

% diff mat
W = zeros(n);
for j=1:n
    for k=(1):n
        W(j,k) = 1/(xj(j) - xj(k));
    end
end
W(eye(n)==1) = 1;
omegaj = prod(W);

% the matrix 1/(x-x_j)
I = eye(n);
diff = x(:)-xj(:).';
A = 1./diff;

% too close
if min(abs(diff(:)))<thd
    [row,col] = find(abs(diff)<thd);
    for i=1:numel(row)
        A(row(i),:) = I(col(i),:);
    end
end

B = A*diag(omegaj);
w = B./sum(B,2);

vals = w*yj(:);

end