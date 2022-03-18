function [coeffs,dd_coeffs] = solve_prony_for_coeffs(x,ord_v,moments)
	n_coeffs = sum(ord_v);
	V = cvand(x,ord_v,n_coeffs);
	coeffs = V\moments(1:n_coeffs);
    dd_coeffs = get_dd_coeffs(x,ord_v,coeffs);
end

