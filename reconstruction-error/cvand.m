function V = cvand(x,ord_v,nrows)
  % Returns confluent Vandermonde matrix
  n_coeffs = sum(ord_v);
  np = length(x);
  V = zeros(n_coeffs,n_coeffs);

  for i=0:nrows-1
      cur_col = 1;
      for j=1:np
	  for k=0:ord_v(j)-1
	      if (k <= i)
		  V(i+1,cur_col) = ff(i,k)*x(j)^(i-k);
	      end
	      cur_col = cur_col+1;
	  end
      end
  end
end