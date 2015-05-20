function [ Nout ] = coxder( N, knotSpan, k)
% COXDER computes derivatives of basis splines by Cox-de Boor recursion. It
% is assumed that N is the output of COX(x,pk), i.e., its first row
% contains the function values at x. This function fills out the second row
% with the first derivative. The 2nd to pth derivatives are not yet
% implemented.
% N: (p+1)x(p+1) matrix, return value of COX(x,p,k).
% knotSpan: knot span of the abscissae, return value of COX(x,p,k)
% k: knot vector
% Nout: overwritten array N

Nout = N;
p = size(N,1)-1;

% get supporting knots 'subKnot'
subKnot = zeros(1,2*p+1);
for s=1:(2*p+1)    
    subKnot(s) = k(knotSpan-p+s-1);
end

s = 2;
    
for t=1:p
    
    if((subKnot(t+p)-subKnot(t))~=0)   
        A = p/(subKnot(t+p)-subKnot(t));
    else
        A = 0;
    end
    
    if((subKnot(t+p+1)-subKnot(t+1))~=0)
        B = p/(subKnot(t+p+1)-subKnot(t+1));
    else
        B = 0;
    end
    
    Nout(s,t) = A*N(s,t)-B*N(s,t+1);
    
end

if((subKnot(p+1+p)-subKnot(p+1))~=0)
    A = p/(subKnot(p+1+p)-subKnot(p+1));
else
    A = 0;
end

Nout(s,p+1) = A*N(s,p+1);


