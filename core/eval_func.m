function [ g, gu, gv ] = eval_func( u, v, ku, kv, p, cps )
% EVAL_FUNC evaluates a function on a B-spline patch defined by its control
% values.
% u,v: abscissae
% ku, kv: knot vector
% p: polynomial degrees
% cps: control values of the function

[Nu,su] = cox(u,p(1),ku);
[Nv,sv] = cox(v,p(2),kv);
Nu = coxder(Nu,su,ku);
Nv = coxder(Nv,sv,kv);

g = 0;
gu = 0;
gv = 0;

for i=1:p(1)+1
    
    for j=1:p(2)+1
      
        I = su - p(1) - 1 + i;          
        J = sv - p(2) - 1 + j;
        
        g = g + cps(I,J)*Nu(1,i)*Nv(1,j);
        gu = gu + cps(I,J)*Nu(2,i)*Nv(1,j);
        gv = gv + cps(I,J)*Nu(1,i)*Nv(2,j);
        
    end
 
end

end

