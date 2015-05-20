function [ h ] = plot_gradient_field( gu, gv )
% PLOT_GRADIENT_FIELD visualizes a gradient field in false coloring.
% gu, gv: gradient field
% h: handle to the figure created

mask = ones(size(gu));
mask(find(gu.^2+gv.^2==0)) = 0;

gut = gu(:);
gvt = gv(:);

nx = -gut;
ny = -gvt;
nz = ones(size(gut));

nnorm = sqrt(nx.^2+ny.^2+nz.^2);
nx = nx./nnorm;
ny = ny./nnorm;
nz = nz./nnorm;

C = zeros([size(gu),3]);
C(:,:,1) = mask.*reshape(nx,size(gu));
C(:,:,2) = mask.*reshape(ny,size(gu));
C(:,:,3) = mask.*reshape(nz,size(gu));

z = zeros(size(gu));

h=figure('Name','Gradient field');
set(h,'Color','black');
surf(z,C,'EdgeColor','none','FaceColor','interp');
view(0,90);
axis equal;
axis tight;

end

