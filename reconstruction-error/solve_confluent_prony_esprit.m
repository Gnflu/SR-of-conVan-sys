function  [reconstructed_jumps, reconstructed_coeffs] = solve_confluent_prony_esprit(np,ord,measurements)
	reconstructed_jumps = estimate_jumps(measurements, np,ord);
	multiplicity_structure = (ord)*ones(1,np);

    [reconstructed_coeffs,dd_coeffs] = ...
        solve_prony_for_coeffs(reconstructed_jumps,multiplicity_structure,measurements);

end


function estimated_jumps = estimate_jumps(data_vector, njumps, ord)

	% Build the data matrix
	L = njumps*(ord)+1;
	N = 2*L+1;
	dm = data_matrix(data_vector,N,L);
	[u,w,v]=svd(dm,0);
%	dmr = u*mpower(w(1:L,1:L),1/2);
%	f = dmr(1:L-1,:) \ dmr(2:L,:);
	dmr = u(1:L,1:L);
	f = dmr(1:L-1,1:L-1) \ dmr(2:L,1:L-1);
	
	roots = eig(full(f));
	roots = roots(1:L-1);

    roots_real = [real(roots), imag(roots)];

	idx = kmeans(roots_real,njumps);
	estimated_jumps = zeros(njumps,1);
	for j=1:njumps
		estimated_jumps(j) = mean(roots(idx==j));
	end
end

function dm = data_matrix(seq, N,L)
	if (length(seq) < N)
		error ('Error: there should be at least %d coefficients',N);
	end
    d1 = hankel(seq(1:L),seq(L:N));
    %d2 = conj(rot90(d1,2));
    %dm = [d1 d2];

	dm = d1;
end