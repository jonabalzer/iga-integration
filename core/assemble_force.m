function [ force ] = assemble_force( ku, kv, p, gu, gv )
% ASSEMBLE_FORCE assembles the right-hand side of the discrete Poisson 
% equation characterizing the minimum of the gradient error functional.
% ku, kv: knot vectors 
% p: polynomial degrees
% gu, gv: gradient field to integrate
% force: force vector

disp('Assembling force vector...');

% get finite elements
eu = get_element_vector(ku,p(1));
ev = get_element_vector(kv,p(2));

% get some stats
nb = (p(1)+1)*(p(2)+1);                     % no of local basis functions
n = [ncp(ku,p(1)),ncp(kv,p(2))];            % no of control points per dim
ngps = [eu(2),ev(2)]-0.5;                   % no of Gauss points

% allocate temp and output arrays
force = zeros(n(1)*n(2),1);

% iterate through the elements
for e=1:length(eu)-1
    
    for f=1:length(ev)-1
        
        % evaluate at Gauss points = image pixels
        for gpi=1:ngps(1)
            
            for gpj=1:ngps(2)
                
                [Nu,su] = cox(eu(e)-0.5+gpi,p(1),ku);
                [Nv,sv] = cox(ev(f)-0.5+gpj,p(2),kv);
                Nu = coxder(Nu,su,ku);
                Nv = coxder(Nv,sv,kv);
                
                for i=1:(p(1)+1)
                    
                    for j=1:(p(2)+1)
                        
                        gm = [gu(eu(e)-0.5+gpi,ev(f)-0.5+gpj);gv(eu(e)-0.5+gpi,ev(f)-0.5+gpj)];
                
                        cube(i,j,(gpi-1)*ngps(2)+gpj) = gm(1)*Nu(2,i)*Nv(1,j) + gm(2)*Nu(1,i)*Nv(2,j);
                        
                    end
                    
                end
                
            end
            
        end
        
        % actual integration
        for i=1:(p(1)+1)
                    
           for j=1:(p(2)+1)
                
                int = 0;
                
                for gpi=1:ngps(1)
                    
                    for gpj=1:ngps(2)
                        
                        int = int + cube(i,j,(gpi-1)*ngps(2)+gpj);
                        
                    end
                    
                end
                
                I = e - 1 + i;
                J = f - 1 + j;
                A = (I-1)*n(2) + J;
                
                force(A) = force(A) + int;
                
            end
            
        end
        
        
    end
    
    
end

force = [force;0];

end
