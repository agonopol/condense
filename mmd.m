function [d]= mmd(X, Y)
    n = size(X,1);
    m = size(Y, 1);
    s = sigma([X;Y]);
    
    Kxx = grbf(X, X, s);
    Kyy = grbf(Y, Y, s);
    Kxy = grbf(X, Y, s);
    d = (1.0 / (n ^ 2.0)) *  sum(sum(Kxx)) ... 
             - (2.0 / ( n * m )) * sum(sum(Kyy)) ... 
             + (1.0 / ( m ^ 2.0)) * sum(sum(Kxy));
%     d = sum(sum(Kxx)) + sum(sum(Kyy)) - (sum(diag(Kxx)) + sum(diag(Kyy))) ...
%         - 2 * (sum(sum(Kxy)) - (sum(diag(Kxy))));
%     d = d / (n*(m-1));
    d = sqrt(abs(d));
end

function d = grbf(p, q, sigma)
    n = size(p, 1);
    m = size(q, 1);
    
    k1 = sum(p .* p, 2);
    r = repmat(k1', m, 1)';
    
    k2 = sum(q .* q, 2);
    s = repmat(k2', n, 1);
    
    h = r + s;
    h = h - 2 * p * q';
    d = exp(-1*h/(2*sigma^2.0));
end

function [s] = sigma(M)
 D = squareform(pdist(M));
 s = mean2(D) - std2(D);  
end


