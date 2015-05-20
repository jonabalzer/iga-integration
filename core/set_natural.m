function [ L,f ] = set_natural( L, f, gu, gv, ku, kv, p, cpx, cpy, cpz )
% SET_NATURAL enforces the natural boundary condition along the entire
% boundary of the spline patch.
% L: stiffness matrix
% f: force vector
% gu, gv: gradient data
% ku, kv: knot vectors
% p: polynomial degrees
% cpx, cpy, cpz: control points

disp('Setting boundary conditions...');

n = [ncp(ku,p(1)),ncp(kv,p(2))];

[u,v] = meshgrid(1:size(gu,1),1:size(gu,2));
guint = interp2(u,v,gu,cpx,cpy,'*linear',0)';
gvint = interp2(u,v,gv,cpx,cpy,'*linear',0)';

for I=1:size(cpx,1)
    
    for J=1:size(cpx,2)

        fval = 0;
        dir = [0;0];
 
        if(I==1 && J==1)
            dir = [-1;-1]/sqrt(2);
            fval = dir(1)*guint(I+1,J+1) + dir(2)*gvint(I+1,J+1);
        end
        
        if(I==1 && J>1 && J<size(cpx,2))
            dir = [-1;0];
            fval = dir(1)*guint(I+1,J) + dir(2)*gvint(I+1,J);
        end
        
        if(I==1 && J==size(cpx,2))
            dir = 0.5*[-1;1]/sqrt(2);
            fval = dir(1)*guint(I+1,J-1) + dir(2)*gvint(I+1,J-1);
        end
        
        if(I==size(cpx,1) && J==1)
            dir = 0.5*[1;-1]/sqrt(2);
            fval = dir(1)*guint(I-1,J+1) + dir(2)*gvint(I-1,J+1);
        end
        
        if(I==size(cpx,1) && J>1 && J<size(cpx,2))
            dir = [1;0];
            fval = dir(1)*guint(I-1,J) + dir(2)*gvint(I-1,J);
        end
        
        if(I==size(cpx,1) && J==size(cpx,2))
            dir = 0.5*[1;1]/sqrt(2);
            fval = dir(1)*guint(I-1,J-1) + dir(2)*gvint(I-1,J-1);
        end
        
        if(J==1 && I>1 && I<size(cpx,1))
            dir = [0;-1];
            fval = dir(1)*guint(I,J+1) + dir(2)*gvint(I,J+1);
        end
        
        if(J==size(cpx,2) && I>1 && I<size(cpx,1))
            dir = [0;1];
            fval = dir(1)*guint(I,J-1) + dir(2)*gvint(I,J-1);
        end
        
     
        if(norm(dir,2)~=0)
        
            cp = [cpx(I,J),cpy(I,J),cpz(I,J)];
            cpi = [I,J];
            
            L = set_neumann_bc(L,cpi,cp,-dir,ku,kv,p);
            
            f((I-1)*n(2) + J) = -0.5*fval;
                        
        end

    end
    
end




end

