%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Copyright (C) 2012  Jonathan Balzer
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init
clear all;
close all;
clc;
path(path, './core');

% load data
load('./data/sine.mat');
npx = size(gu);

% set spline parameters
p = [2,2];
scales = [4,4];

% add noise
sigma = 0;
nu = normrnd(0,sigma,npx);
nv = normrnd(0,sigma,npx);
gum = gu + nu;
gvm = gv + nu;

% corrupt data
width = 0;
density = 0.001; 
mask = imnoise(zeros(size(gu)),'salt & pepper',density);
h = ones(width,width);
mask = imfilter(mask,h);
mask(find(mask>1)) = 1;
mask = ones(size(mask)) - mask;
gum = gum.*mask;
gvm = gvm.*mask;

% initialize the B-spline patch
ku = init_knot_vector(npx(1),scales(1),p(1));   % knot vectors
kv = init_knot_vector(npx(2),scales(2),p(2));
[cpx,cpy,cpz] = init_cps(ku,kv,p);              % control point grid

% assemble Poisson equation
t = cputime;
L = assemble_lb(ku,kv,p);
f = assemble_force(ku,kv,p,gum,gvm);
disp(['Assembly finishes after ',num2str(cputime-t),' s.']);

% solve the system
disp(['Solving sparse linear system of size ',int2str(size(L,1)-1),'x',int2str(size(L,2)-1),'...']);
t = cputime;
sol = L\f;                                                         % direct solver
% sol = pcg(L,f,1e-6,500);                                         % optional: pre-conditioned iterative CG
disp(['Solved after ',num2str(cputime-t),' s.']);

% remove normalization = mean-value condition
sol = sol(1:length(f)-1);                                          
cpz = reshape(sol,size(cpx))';                                     

% evaluate function and gradients at grid points
[u,v] = ndgrid(1:1:npx(1),1:1:npx(2));
gr = zeros(size(u));
gur = gr;
gvr = gr;

for i=1:size(u,1)
     
     for j=1:size(u,2)
         
         [gr(i,j),gur(i,j),gvr(i,j)] = eval_func(u(i,j),v(i,j),ku,kv,p,cpz);
       
     end
     
end

% L2-error
locerror = sqrt((gu-gur).*(gu-gur)+(gv-gvr).*(gv-gvr));
error = norm((gu-gur)+(gv-gvr),2)/sqrt(npx(1)*npx(2));
disp(['L2-error:',num2str(error)]);

% do some visualizations
plot_gradient_field(gum,gvm);
 
h=figure('Name','Error');
set(h,'Color','black');
hold on
surf(u,v,gr,locerror,'EdgeColor','none','FaceColor','interp');
caxis([min(min(locerror)),max(max(locerror))]);
cbar = colorbar('Location','West');
cpos = get(cbar,'Position');
cpos(3)=1.5*cpos(3);
cpos(4)=0.25*cpos(4);
cpos(2)=cpos(2)+0.12;
cpos(1)=cpos(1)-0.1;
set(cbar,'Position',cpos);
axis equal;
view(-40,21); 
axis tight
camlight; 
lighting phong;
set(gca,'Visible','off')
camzoom(1.75);
camdolly(0.1,-0.2,0)

h=figure('Name','Spline reconstruction');
set(h,'Color','black');
hold on;
lighting phong;
surf(u,v,gr,'EdgeColor','none','FaceColor',[0.7,0.7,0.7]);
mesh(cpx,cpy,cpz,'FaceColor','none','EdgeColor',[1,0,0],'Marker','o','Linewidth',2,'EdgeLighting','none','MarkerFaceColor',[1,0,0],'MarkerSize',14);
axis equal;
axis tight
set(gca,'Visible','off')
view(-40,21);
lighting phong;
camlight; 
camzoom(1.7);
camdolly(0.1,0.2,0)
 
% if ground truth available, plot it
if(norm(g,2)>0)
    
    h=figure('Name','Ground truth surface');
    set(h,'Color','black');
    hold on
    surf(u,v,g,'EdgeColor','none','FaceColor',[0.7,0.7,0.7]);
    axis equal;
    view(-40,21);
    axis tight
    camlight;
    lighting phong;
    set(gca,'Visible','off')
    camzoom(1.75);
    camdolly(0.1,0.2,0)
    
end


