function [ n ] = ncp( k, p )
% NCP calculates the number of control points in one dimension. 
% k: knot vector
% p: polynomial degree

n = length(k)-p-1;

end

