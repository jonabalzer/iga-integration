function [ e ] = get_element_vector( k, p )
% GET_ELEMENT_VECTOR extracts isogeometric finite elements from a knot
% vector.
% k: knot vector
% p: polynomial degree

e = k(p+1:length(k)-p);

end

