% TODO make better generally
% TODO error handling
% Solve Ax=b using Gaussian elimination without pivoting
function x = naive_gaussian(A, b)
  dim = length(A);
  x = zeros(dim,1);
  
  % forward elimination
  for i = 2:dim
    m = -A(i:dim,i-1) / A(i-1,i-1);
    A(i:dim,:) += m * A(i-1,:);
    b(i:dim,:) += m * b(i-1,:);
  endfor
  
  % backward substitution
  x(dim,:) = b(dim,:) / A(dim,dim);
  for i = dim:-1:1
    x(i,:) = b(i,:) - A(i,i+1:dim) * x(i+1:dim,:);
    x(i,:) /= A(i,i);
  endfor
endfunction
