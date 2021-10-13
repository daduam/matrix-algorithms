% TODO make better generally
% TODO error handling
% Solve Ax=b using Doolittle LU decomposition
function x = doolittle_lu(A,b)
  dim = length(A);
  L = eye(dim);
  U = A;
  
  % building L and U matrices
  for j = 1:dim-1
    for i = j+1:dim
      L(i,j) = U(i,j) / U(j,j);
      U(i,j:dim) -= L(i,j) * U(j,j:dim);
    endfor
  endfor
  
  % LY = b then Ux = Y
  % LY = b. find Y using forward substitution
  Y = zeros(dim, 1);
  Y(1) = b(1) / L(1,1);
  for i = 2:dim
    Y(i) = b(i) - L(i,1:i-1) * Y(1:i-1);
    Y(i) /= L(i,i);
  endfor
  
  % Ux = Y. find x using back substitution
  x = zeros(dim, 1);
  x(dim,:) = Y(dim,:) / U(dim,dim);
  for i = dim:-1:1
    x(i,:) = Y(i,:) - U(i,i+1:dim) * x(i+1:dim,:);
    x(i,:) /= U(i,i);
  endfor
endfunction
