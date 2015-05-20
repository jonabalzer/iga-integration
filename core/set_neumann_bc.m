function [ L ] = set_neumann_bc( L, cpi, cp, dir, ku, kv, p )
% SET_NEUMANN_BC modifies the stiffness matrix as to enforce the natural or 
% free boundary condition at a given set of points. The desired value has to
% be set at an appropriate location in the force vector.
% L: stiffness matrix
% cpi: control points at whose Greville abscissae the condition is imposed
% cp: control points
% dir: direction of derivative
% ku, kv: knot vectors
% p: polynomial degree

% get some stats
n = [ncp(ku,p(1)),ncp(kv,p(2))];            % no of control points per dim

% convert to lexicographic index
A = (cpi(1)-1)*n(2) + cpi(2);

% erase original 
L(A,:) = 0;

% evaluate at Greville abscissae
[Nu,su] = cox(cp(1),p(1),ku);
[Nv,sv] = cox(cp(2),p(2),kv);
Nu = coxder(Nu,su,ku);
Nv = coxder(Nv,sv,kv);

for i=1:p(1)+1
    
    for j=1:p(2)+1
      
        I = su - p(1) - 1 + i;          
        J = sv - p(2) - 1 + j;
        B = (I-1)*n(2) + J;
        
        L(A,B) = dir(1)*Nu(2,i)*Nv(1,j) + dir(2)*Nu(1,i)*Nv(2,j);
               
    end
 
end


end

