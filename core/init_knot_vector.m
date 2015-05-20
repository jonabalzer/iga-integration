function [ k ] = init_knot_vector( nopx, scale, p )
% INIT_KNOT_VECTOR creates a knot vector whose entries are distributed
% equally over a Cartisian image grid.
% nopx: number of pixels
% scale: density of knots
% p: polynomial degree
% k: knot vector

ld = log(nopx)/log(2);

if(mod(ld,1)~=0)
    error('Image size not a power of 2!');  
end

if(scale>ld || scale <0)
    error(['Scale must be between 0 and ',int2str(ld),'!']);  
end

noElements = 2^(ld-scale);

k = [zeros(1,p),0:noElements:nopx,nopx*ones(1,p)]+0.5;

end

