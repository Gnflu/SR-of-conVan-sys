function res = ff(x,y)
% Calculate falling factorial x^{(y)}
    if (x==0)
        res = 1;
    else
        res = prod(x:-1:max(x-y+1,1));
    end
end