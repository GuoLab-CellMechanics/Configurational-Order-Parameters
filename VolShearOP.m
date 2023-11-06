%
%
% Example code of volumetric and shear order parameters for multicellular systems
% See our paper for more details: https://www.pnas.org/doi/abs/10.1073/pnas.2109168118
% 
% 

clear; clf;

% load the files
load('coordinates.mat');
image=imread('Composite2.jpg');

% Delaunay triangulation
po=[x;y];
ind = delaunayn(po');

% get the referential area
cell_num      = length(x);
total_area    = (max(x)-min(x))*(max(y)-min(y));
area_baseline = total_area/2/cell_num;

for i=1:size(ind,1)
    
    % get the index of the three vertices of a triangle
    j=ind(i,1);
    p=ind(i,2);
    q=ind(i,3);
    
    % the triangle matrix
    dx=[x(p)-x(j), x(q)-x(j);
        y(p)-y(j), y(q)-y(j);];
    
    % the unit triangle matrix for reference
    dX=[1,1/2;
        0,sqrt(3)/2;];
    
    % calculate deformation
    temp_F=dx/dX;               % deformation gradient
    temp_E=temp_F'*temp_F;      % Cauchy-Green deformation tensor
    I1=trace(temp_E);           % The first invariant
    jacobi=abs(det(temp_F));    % Jacobian
    
    % calculate the local volume and shear order
    vol(i)=log(1/2*jacobi/area_baseline*sqrt(3)/2);
    shear(i)=log(I1/jacobi-2);
end

% get coordinates of each triangle for plots
for e=1:size(ind,1)
    P_x(:,e)=x(ind(e,:));
    P_y(:,e)=y(ind(e,:));
end

% plot figures
figure(1)
imshow(image)
hold on
patch(P_x,P_y,vol,'facecolor','none','edgecolor',[255 0 0]/255,'edgealpha',0.5,'linewidth',3)
set(gca,'YDir','normal')
xlim([256 768])
ylim([256 768])

figure(2)
patch(P_x,P_y,vol,'EdgeColor','none','edgecolor',[255 0 0]/255,'edgealpha',0.2)
cb=colorbar;
caxis([-2,2])
axis equal
axis off
title('lg(vol)')
xlim([256 768])
ylim([256 768])

figure(3)
hold on;
patch(P_x,P_y,shear,'EdgeColor','none','edgecolor',[255 0 0]/255,'edgealpha',0.2)
cb=colorbar;
caxis([-2,2])
axis equal
axis off
title('lg(shear)')
xlim([256 768])
ylim([256 768])

% calculate the global order parameters
phi_vol   = mean(vol.^2,'omitnan');
phi_shear = mean(shear,'omitnan');
sprintf('%s%03f%s%03f','phi_vol = ',phi_vol,',phi_shear = ',phi_shear)

