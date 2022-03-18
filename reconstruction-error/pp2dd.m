function A = pp2dd(x,ord_v)
  n_coeffs = sum(ord_v);
  np = length(x);
  VV = cvand(x,ord_v,n_coeffs)';
  A = zeros(n_coeffs,n_coeffs);
  curr_col = 1;
  for pp=1:np
      for dd=0:ord_v(pp)-1
	  % Fill the column
	  idxs = 1:curr_col;
	  U = VV(idxs,idxs);
	  du = det(U);
	  UU = U(:,1:end-1);
	  for i=1:curr_col
	      remaining_rows = setdiff(idxs,i);
	      if (~isempty(remaining_rows))
		  UUU = UU(remaining_rows,:);
		  duu = det(UUU);
	      else
		  duu = 1;
	      end
	      A(i,curr_col) = (1-2*mod(i+curr_col,2))*duu/du;
	  end
	  curr_col = curr_col+1;
      end
  end
end