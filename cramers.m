% TODO make better generally
% TODO error handling
% Solve Ax=b using Cramer's Rule
function x = cramers(A, b)
  dim = length(A);
  det_a = det(A);
  x = zeros(dim, 1);
  for j = 1:dim
    Aj = A;
    Aj(:,j) = b;
    x(j) = det(Aj) / det_a;
  endfor
endfunction
