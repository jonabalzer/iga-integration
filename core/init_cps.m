function [ cpx, cpy, cpz ] = init_cps( ku, kv, p )
% INIT_CPS distributes the control points of a B-spline patch equally over
% its domain.
% ku, kv: knot vectors
% p: polynomial degrees
% cpx, cpy, cpz: control points

nu = ncp(ku,p(1));
nv = ncp(kv,p(2));

kumin = min(ku);
kumax = max(ku);
kvmin = min(kv);
kvmax = max(kv);

hu = (kumax-kumin)/(nu-1);
hv = (kvmax-kvmin)/(nv-1);

[cpx,cpy] = ndgrid(kumin:hu:kumax,kvmin:hv:kvmax);

cpz = zeros(size(cpx));

end

