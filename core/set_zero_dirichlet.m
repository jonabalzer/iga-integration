function [ L, f ] = set_zero_dirichlet( L, f, ku, kv, p, cpx, cpy, cpz  )
% SET_ZERO_DIRICHLET alters the linear system discretizing the Poisson 
% equation such that the boundary of the solution is fixed.
% L: stiffness matrix
% f: force vector
% ku, kv: knot vectors
% p: polynomial degrees
% cpx, cpy, cpz: control points

disp('Setting boundary conditions...');

% remove normalization
ls = size(L);
L(ls(1),:) = 0;
L(:,ls(2)) = 0;
L(ls(1),ls(2)) = 1;

f(length(f)) = 0;

n = [ncp(ku,p(1)),ncp(kv,p(2))];

for I=1:size(cpx,1)
       
    for J=1:size(cpx,2)
        
        if(I==1 || J==1 || I==size(cpx,1) || J==size(cpx,2))
                        
            f((I-1)*n(2) + J) = 0;
            cp = [cpx(I,J),cpy(I,J),cpz(I,J)];
            cpi = [I,J];
            
           L = set_dirichlet_bc(L,cpi,cp,ku,kv,p);
                     
        end
       
    end
    
end

end

