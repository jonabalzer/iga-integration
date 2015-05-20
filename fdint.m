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
[u,v] = ndgrid(1:1:npx(1),1:1:npx(2));

% add noise
sigma = 1.0;
nu = normrnd(0,sigma,npx);
nv = normrnd(0,sigma,npx);
gum = gu + nu;
gvm = gv + nu;

% add clutter
width = 20;
density = 0.001; 
mask = imnoise(zeros(size(gu)),'salt & pepper',density);
h = ones(width,width);
mask = imfilter(mask,h);
mask(find(mask>1)) = 1;
mask = ones(size(mask)) - mask;
gum = gum.*mask;
gvm = gvm.*mask;

% assemble Poisson equation
t = cputime;
[Ld,f] = fd_discretization(gum,gvm);
disp(['Assembly finishes after ',num2str(cputime-t),' s.']);

% solve sparse system
disp(['Solving sparse linear system of size ',int2str(size(Ld,1)),'x',int2str(size(Ld,2)),'...']);
t = cputime;
gr = Ld\f;
disp(['Solved after ',num2str(cputime-t),' s.']);
gr = gr(1:length(f)-1);
gr = reshape(gr,size(gu));

% L2 error
[gvr,gur] = gradient(g);
locerror = sqrt((gu-gur).*(gu-gur)+(gv-gvr).*(gv-gvr));
error = norm((gu-gur)+(gv-gvr),2)/sqrt(npx(1)*npx(2));
disp(['L2-error:',num2str(error)]);

% visualization
plot_gradient_field(gum,gvm);

h=figure('Name','Error');
set(h,'Color','black');
hold on
surf(u,v,gr,locerror,'EdgeColor','none');
axis equal;
view(-40,21); 
axis tight
camlight; 
lighting phong;
set(gca,'Visible','off')
camzoom(1.75);
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
