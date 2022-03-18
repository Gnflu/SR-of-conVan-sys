clc
clear all
close all
format long

s = 7; % # of nodes
l = 3; % # of clustered nodes
nt = 1000; % # of test
hs = 0.0005:0.00001:0.001; % size of cluster
Ns = 100:1:200; % # of samples

plot_single_cluster_sv_rate(s,l,hs,Ns,nt);

% create a single clustered configuration (normally distributed)
function x = clustered_config(p,d,h)
x = [];
th = h/(p-1); %tau*h
    for j=1:1:p
        x(j) = (j-1)*th;
    end
    for j=(p+1):1:d
        x(j) = (p-1)*th + ((j-p)*((pi - (p-1)*th))/(d-p+1));
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
    
function plot_single_cluster_sv_rate(s,l,hs,Ns,nt)
    SRFs = ones([1 nt]);
    sv = ones([1 nt]);
    for i = 1:nt
        randomIndex = randi(length(hs), 1);
        h = hs(randomIndex);
        randomIndex = randi(length(Ns), 1);
        N = Ns(randomIndex);
        x = clustered_config(l,s,h);
        U1 = create_confluent_van(x, N);
        U = U1./sqrt(N);
        S = svd(U);
        sv(i) = S(end);
        SRFs(i) = (N*(h/(s-1)))^(-1);
    end
    %Plot figure
    p = polyfit(log10(SRFs),log10(sv),1);
    y = polyval(p,log10(SRFs));
    figure
    abb=scatter(log10(SRFs),log10(sv),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7]);
    hold on
    ab=plot(log10(SRFs),log10((SRFs).^(1-(2*l)))-2,'LineWidth',2);
    grid on
    titleInfo = sprintf('Single cluster l = %d, s = %d', l, s);
    legend('$\sigma_{min}$','$SRF^{1-2{\ell}}$','Interpreter','latex');
    title(titleInfo,'Interpreter','latex')
    xlabel('SRF','Interpreter','latex')
end


