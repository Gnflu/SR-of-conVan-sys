clc
clear all
close all
format long

s = 7; % # of nodes
l1 = 2; % # of nodes in cluster 1
l2 = 3; % # of nodes in cluster 2
nt = 1000; % # of test
hs = 0.0005:0.00001:0.001; % size of cluster
Ns = 100:1:200; % # of samples


plot_multi_cluster_sv_rate(s,l1,l2,hs,Ns,nt);


% create a multi clustered configuration (normally distributed)
function x = two_clustered_config(p1,p2,d,h1,h2)
x = [];
th1 = h1/(p1-1); %tau*h
th2 = h2/(p2-1); %tau*h
    for j=1:1:p1
        x(j) = (j-1)*th1;
    end
    for j=(p1+1):1:(d-p2+1)
        x(j) = (p1-1)*th1 + ((j-p1)*((pi - (p1-1)*th1 - (p2-1)*th2))/(d-p1-p2+2));
    end
    for j=(d-p2+2):1:d
        x(j) = x(d-p2+1) + (d-j+1)*th2;
    end
    x = x - (pi/2);
end



function U = create_confluent_van(clfun, M)
    d = length(clfun);
    G = diag(0:M-1);
    I1 = eye(M);
    L = [I1 G];
    z = exp(1i*clfun);
    V = vanderm(z,M);
    V2 = blkdiag(V,V);
    I2 = eye(d);
    D = diag(z.^-1);
    R = blkdiag(I2,D);
    U = L*V2*R;
end


function V = vanderm(v,N)
    %  V = vanderm(v,N) returns the n*N Vandermonde matrix whose columns are
    %  powers of the vector v, i.e. A(i,j) = v(i)^(j-1), where n is the length
    %  of v. If N is not specified, square matrix is assumed, i.e. N=n.
    if all(size(v)>1)
        error('Input v should be a vector.')
    end
    narginchk(1,2)
    nargoutchk(0,2)
    v = v(:); % make v a column
    n = length(v);
    if nargin < 2
        N = n; 
    end
    X = [ones(n,1) repmat(v,[1 N-1])];
    A = cumprod(X,2);
    V = A.';
end 
    

function plot_multi_cluster_sv_rate(s,l1,l2,hs,Ns,nt)
    for i = 1:nt
        randomIndex1 = randi(length(hs), 1);
        randomIndex2 = randi(length(hs), 1);
        h1 = hs(randomIndex1);
        h2 = hs(randomIndex2);
        randomIndex = randi(length(Ns), 1);
        N = Ns(randomIndex);
        x = two_clustered_config(l1,l2,s,h1,h2);
        U1 = create_confluent_van(x, N);
        U = U1./sqrt(N);
        S = svd(U);
        sv(i) = S(end);
   
        h = min(h1,h2);
        l = max(l1,l2);
        SRFs(i) = ((N)*h/(l-1))^(-1);
    end
    %Plot figure
    p = polyfit(log10(SRFs),log10(sv),1);
    y = polyval(p,log10(SRFs));
    figure
    abb=scatter(log10(SRFs),log10(sv),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7]);
    hold on
    ab=plot(log10(SRFs),log10((SRFs).^(1-(2*l)))-5.2,'LineWidth',2);
    grid on
    legend('$\sigma_{min}$','$SRF^{1-2{\ell}_{max}}$','Interpreter','latex');
    titleInfo = sprintf('Multi cluster l1 = %d, l2 = %d, s = %d', l1, l2, s);
    title(titleInfo,'Interpreter','latex')
    xlabel('SRF','Interpreter','latex')
end
  
