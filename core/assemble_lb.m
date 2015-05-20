function [ L ] = assemble_lb( ku, kv, p )
% ASSEMBLE_LB assembles the stiffness matrix for the Poisson problem, i.e.,
% the discrete Laplace-Beltrami operator. 
% ku,kv: knot vectors
% p: polynomial degrees

disp('Assembling stiffness matrix...');

% get finite elements
eu = get_element_vector(ku,p(1));
ev = get_element_vector(kv,p(2));

% get some stats
nb = (p(1)+1)*(p(2)+1);                     % no of local basis functions
n = [ncp(ku,p(1)),ncp(kv,p(2))];            % no of control points per dim
ngps = [eu(2),ev(2)]-0.5;                   % no of Gauss points

% allocate temp and output arrays
rowInds = zeros(1,(length(eu)-1)*(length(ev)-1)*nb*nb); %+n(1)*n(2)+1);    
colInds = zeros(1,(length(eu)-1)*(length(ev)-1)*nb*nb); %+n(1)*n(2)+1);
values = zeros(1,(length(eu)-1)*(length(ev)-1)*nb*nb); %+n(1)*n(2)+1);
cube = zeros(nb,nb,ngps(1)*ngps(2));
Nau = zeros(1,nb);
Nav = zeros(1,nb);

counter = 1;  % matrix entry counter

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
                        
                        Nau((i-1)*(p(2)+1)+j) = Nu(2,i)*Nv(1,j);
                        Nav((i-1)*(p(2)+1)+j) = Nu(1,i)*Nv(2,j);
                        
                    end
                    
                end
                
                for a=1:nb
                    
                    for b=1:nb
                        
                        cube(a,b,(gpi-1)*ngps(2)+gpj) = Nau(a)*Nau(b) + Nav(a)*Nav(b);
                        
                    end
                    
                    
                end
                
            end
            
        end
        
        % perform the actual integration
        for a=1:nb
            
            for b=1:nb
                
                int = 0;
                
                for gpi=1:ngps(1)
                    
                    for gpj=1:ngps(2)
                        
                        int = int + cube(a,b,(gpi-1)*ngps(2)+gpj);
                        
                    end
                    
                end
                
                Ia = e + floor((a-1)/(p(2)+1)); % e - 1 + i
                Ja = f + mod(a-1,p(2)+1);
                A = (Ia-1)*n(2) + Ja;
                
                Ib = e + floor((b-1)/(p(2)+1));
                Jb = f + mod(b-1,p(2)+1);
                B = (Ib-1)*n(2) + Jb;
                
                rowInds(counter) = A;
                colInds(counter) = B;
                values(counter) = int;
                counter = counter + 1;
                
            end
            
        end
        
        
    end
    
end

% initialize sparse matrix
L = sparse(rowInds,colInds,values);
ls = size(L);

% add zero-mean condition
L(ls(1)+1,:) = 1;
L(:,ls(2)+1) = 1;
L(ls(1)+1,ls(2)+1) = 1;


end

