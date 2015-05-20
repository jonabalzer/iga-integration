function [ L, f ] = fd_discretization( gu, gv )
% FD_DISCRETIZATION performs finite-difference approximation of the Euler
% equation under natural boundary conditions of the variational integration
% problem.
% gu,gv: gradient data
% L: discrete Laplacian
% f: discrete divergence of the gradient field

f = zeros(size(gu,1)*size(gu,2),1);

n = size(gu,1);
noEntries = 2*(4*n-2) + 4*(n^2-4*n+2);

rowInds = zeros(1,noEntries);
colInds = zeros(1,noEntries);
val = zeros(1,noEntries);

counter = 1;

for i=2:(size(gu,1)-1)
    
    for j=2:(size(gu,2)-1)
        
        A = (j-1)*size(gu,1) + i;
        
        rowInds(counter) = A;
        colInds(counter) = A;
        val(counter) = -4;
        counter = counter + 1;
        
        rowInds(counter) = A;
        colInds(counter) = (j-2)*size(gu,1) + i;
        val(counter) = 1;
        counter = counter + 1;
        
        rowInds(counter) = A;
        colInds(counter) = j*size(gu,1) + i;
        val(counter) = 1;
        counter = counter + 1;
        
        rowInds(counter) = A;
        colInds(counter) = (j-1)*size(gu,1) + i - 1;
        val(counter) = 1;
        counter = counter + 1;
        
        rowInds(counter) = A;
        colInds(counter) = (j-1)*size(gu,1) + i + 1;
        val(counter) = 1;
        counter = counter + 1;
        
        f(A) = 0.5*(gu(i+1,j)-gu(i-1,j)) + 0.5*(gv(i,j+1)-gv(i,j-1));   % divergence
        
        
    end
    
end


j = 1;
for i=2:size(gu,1)-1
    
    A = (j-1)*size(gu,1) + i;
    
    rowInds(counter) = A;
    colInds(counter) = A;
    val(counter) = 1;
    counter = counter + 1;
    
     f(A) = -gv(i,j);

    rowInds(counter) = A;
    colInds(counter) = j*size(gu,1) + i;
    val(counter) = -1;
    counter = counter + 1;
    
end

j = size(gu,2);
for i=2:size(gu,1)-1
    
    A = (j-1)*size(gu,1) + i;
    
    rowInds(counter) = A;
    colInds(counter) = A;
    val(counter) = 1;
    counter = counter + 1;
    
    f(A)  = gv(i,j);
    
    rowInds(counter) = A;
    colInds(counter) = (j-2)*size(gu,1) + i;
    val(counter) = -1;
    counter = counter + 1;
    
end

i = 1;
for j=2:size(gu,2)-1
    
    A = (j-1)*size(gu,1) + i;
    
    rowInds(counter) = A;
    colInds(counter) = A;
    val(counter) = 1;
    counter = counter + 1;
    
    f(A) = -gu(i,j);
   
    rowInds(counter) = A;
    colInds(counter) = (j-1)*size(gu,1) + i + 1;
    val(counter) = -1;
    counter = counter + 1;
    
end


i = size(gu,1);
for j=2:size(gu,2)-1
    
    A = (j-1)*size(gu,1) + i;
    
    rowInds(counter) = A;
    colInds(counter) = A;
    val(counter) = 1;
    counter = counter + 1;
    
    f(A) = gu(i,j);
    
    rowInds(counter) = A;
    colInds(counter) = (j-1)*size(gu,1) + i - 1;
    val(counter) = -1;
    counter = counter + 1;
    
end

% NW
i=1;
j=1;

A = (j-1)*size(gu,1) + i;

rowInds(counter) = A;
colInds(counter) = A;
val(counter) = 1;
counter = counter + 1;

f(A) = (-gu(i,j)-gv(i,j))/sqrt(2);
      
rowInds(counter) = A;
colInds(counter) = j*size(gu,1) + i + 1;
val(counter) = -1;
counter = counter + 1;

% SW
i=1;
j=size(gu,2);

A = (j-1)*size(gu,1) + i;

rowInds(counter) = A;
colInds(counter) = A;
val(counter) = 1;
counter = counter + 1;

f(A) = (-gu(i,j)+gv(i,j))/sqrt(2);

rowInds(counter) = A;
colInds(counter) = (j-2)*size(gu,1) + i + 1;
val(counter) = -1;
counter = counter + 1;

% NE
i=size(gu,1);
j=1;

A = (j-1)*size(gu,1) + i;

rowInds(counter) = A;
colInds(counter) = A;
val(counter) = 1;
counter = counter + 1;

f(A) = (gu(i,j)-gv(i,j))/sqrt(2);
            
rowInds(counter) = A;
colInds(counter) = j*size(gu,1) + i - 1;
val(counter) = -1;
counter = counter + 1;

% SE
i=size(gu,1);
j=size(gu,2);

A = (j-1)*size(gu,1) + i;

rowInds(counter) = A;
colInds(counter) = A;
val(counter) = 1;
counter = counter + 1;
  
f(A)  = (gu(i,j)+gv(i,j))/sqrt(2);

rowInds(counter) = A;
colInds(counter) = (j-2)*size(gu,1) + i - 1;
val(counter) = -1;
counter = counter + 1;

L = sparse(rowInds,colInds,val);


ls = size(L);
L(ls(1)+1,:)=1;
L(:,ls(2)+1)=1;
L(ls(1)+1,ls(2)+1)=0;

f = [f;0];

end

