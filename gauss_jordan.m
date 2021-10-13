function x = gauss_jordan(A, b)
  dim = length(A);
  
  % augment A with b
  A = [A, b];
  
  for j = 1:dim-1
    for i = 2:dim
      if A(j,j) == 0
        [A(1,:), A(i,:)] = deal(A(i,:), A(1,:));
      endif
    endfor
    % zero out elements below major diagonal
    for i = j+1:dim
      A(i,:) -= A(j,:) * (A(i,j) / A(j,j));
    endfor
  endfor
  % zero out elements above major diagonal
  for j = dim:-1:2
    for i = j-1:-1:1
      A(i,:) -= A(j,:) * (A(i,j) / A(j,j));
    endfor
  endfor
  % make elements on major diagonal unity
  for i = 1:dim
    A(i,:) /= A(i,i);
  endfor

  % result is the last column of [A|b]
  x = A(:,dim+1);
endfunction
