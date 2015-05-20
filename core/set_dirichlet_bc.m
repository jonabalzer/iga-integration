function [ L ] = set_dirichlet_bc( L, cpi, cp, ku, kv, p )
% SET_DIRICHLET_BC modifies the stiffness matrix to enforce the essential 
% boundary condition at a set of points. The desired value has to be set at
% an appropriate location in the force vector.
% L: stiffness matrix
% cpi: control points at whose Greville abscissae the condition is imposed
% cp: control points
% ku, kv: knot vectors
% p: polynomial degree

% get some stats
n = [ncp(ku,p(1)),ncp(kv,p(2))];           

% convert to lexicographic index
A = (cpi(1)-1)*n(2) + cpi(2);

% erase original 
L(A,:) = 0;

% evaluate at Greville abscissae
[Nu,su] = cox(cp(1),p(1),ku);
[Nv,sv] = cox(cp(2),p(2),kv);

for i=1:p(1)+1
    
    for j=1:p(2)+1
      
        I = su - p(1) - 1 + i;          
        J = sv - p(2) - 1 + j;
        B = (I-1)*n(2) + J;
        
        L(A,B) = Nu(1,i)*Nv(1,j);
               
    end
 
end


end

