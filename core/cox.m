function [ N, knotSpan ] = cox( x, p, k )
% COX computes the values of the p+1 basis splines at x. These values are
% put into the first row of a quadratic matrix N, which can be filled up
% with derivative values at x by COXDER(N,knotSpan,k).
% x: abscissa
% p: polynomial degree
% k: knot vector
% N: (p+1)x(p+1) array whose first row contains the result
% knotSpan: knot interval that contains x

N = zeros(p+1,p+1);

knotSpan = -1;

% find knot span
for s=1:(length(k)-1)

    if(x>=k(s) && x<k(s+1))
       
        knotSpan = s;
        break;
        
    elseif (x==k(length(k)))
        
        knotSpan = length(k)-p-1;
        break;
        
    end
    
end

if(knotSpan==-1)
    error('Error! Point out of domain!');
end

N(1,p+1)=1;

% get supporting knots 'subKnot'
subKnot = zeros(1,2*p+1);
for s=1:(2*p+1)    
    subKnot(s) = k(knotSpan-p+s-1);
end

for s=2:(p+1)  %
    
    for t=1:p
        
        % every value N(s,t) depends on N(s-1,t) and N(s-1,t+1) 
        if((subKnot(t+s-1)-subKnot(t))~=0 & (x-subKnot(t))~=0)
        %if((x-subKnot(t))~=0)
            A = (x-subKnot(t))/(subKnot(t+s-1)-subKnot(t));
        else
            A = 0;
        end;
        
        if((subKnot(t+s-1+1)-subKnot(t+1))~=0 & (subKnot(t+(s-1)+1)-x)~=0)
        %if((subKnot(t+(s-1)+1)-x)~=0)     
            B = (subKnot(t+(s-1)+1)-x)/(subKnot(t+s-1+1)-subKnot(t+1));
        else
            B = 0;
        end
            
        N(s,t) = A*N(s-1,t)+B*N(s-1,t+1);
               
    end
    
    if((x-subKnot(p+1))~= 0 & (subKnot(p+1+s-1)-subKnot(p+1))~=0)
    %if((subKnot(p+1+s-1)-subKnot(p+1))~=0)
        A = (x-subKnot(p+1))/(subKnot(p+1+s-1)-subKnot(p+1));
    else
        A = 0;
    end
    
    N(s,p+1) = A*N(s-1,p+1); 
    
end

N = flipud(N);
    
end

