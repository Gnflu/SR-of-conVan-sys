function dd_coeffs = get_dd_coeffs(x,ord_v,coeffs)
  A = pp2dd(x,ord_v);
  dd_coeffs = A\coeffs;
end