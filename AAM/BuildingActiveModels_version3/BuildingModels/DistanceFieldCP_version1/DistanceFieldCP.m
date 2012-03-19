function [Points1n,Points2n]=DistanceFieldCP(Points1,I1,Points2,I2)
% This function DistanceFieldCP can be used to find closest points between
% two datasets. Instead of matching points to points, it matches points to
% interpolated points on the object boundary. Thus the function performs
% 3D point-to-plane or 2D Point-to-contour matching. It only matches points
% with points with approximately the same contour/surface normal.
%
%  [Points1n,Points2n]=DistanceFieldCP(Points1,I1,Points2,I2)
%
% inputs,
%   I1 : A logical image (2D, or 3D) containing a solid version of the object.
%   Points1 : A list of points on the image contour/surface of the object
%   I2 : A logical image (2D, or 3D) containing the object.
%   Points2 : A list of points on the image contour/surface of the object
%
% outputs,
%   Points1n : The (interpolated) closest points of Points1 
%   Points2n : The (interpolated) closest points of Points2
%
%
% Example 2D,
% % Add all the needed folders and functions to Matlab search path
%  functionname='DistanceFieldCP.m';
%  functiondir=which(functionname);
%  functiondir=functiondir(1:end-length(functionname));
%  addpath([functiondir 'functions'])
% 
% % Load two 2D images (logical datatype)
%  I1=imread('images/bat1.png');
%  I2=imread('images/bat2.png');
%
% % Get the contours of the bats
%  [Lines1,Points1]=isocontour(I1,0.5);
%  [Lines2,Points2]=isocontour(I2,0.5);
%
% % Calculate the closest points
%  [Points1n,Points2n]=DistanceFieldCP(Points1,I1,Points2,I2);
%
% % Show the result
%  figure,
%  subplot(1,3,1), imshow(I1);
%  subplot(1,3,2), imshow(I2);
%  subplot(1,3,3), 
%  RGB=zeros([size(I1) 3]);
%  RGB(:,:,1)=double(I1)*0.4; 
%  RGB(:,:,2)=double(I2)*0.4;
%  RGB(:,:,3)=double(I2)*0.4;
%  imshow(RGB); hold on
%  plot(Points1(:,2),Points1(:,1),'b.');
%  plot(Points2(:,2),Points2(:,1),'r.');
%  plot([Points2(:,2) Points2n(:,2)]',[Points2(:,1) Points2n(:,1)]','g');
%  plot([Points1(:,2) Points1n(:,2)]',[Points1(:,1) Points1n(:,1)]','m');
%
% Example 3D,
% % Add all the needed folders and functions to Matlab search path
%  functionname='DistanceFieldCP.m';
%  functiondir=which(functionname);
%  functiondir=functiondir(1:end-length(functionname));
%  addpath([functiondir 'functions'])
%  addpath([functiondir 'polygon2voxel_version1h']);
% % c-code to mex function
%  cd([functiondir 'polygon2voxel_version1h']);
%  mex polygon2voxel_double.c 
%  cd(functiondir);
%
% % Load the 3D jaw triangulated surface-descriptions
%  load('images/jaw3d.mat');
%
% % Make the point lists from the vertices
%  Points1=FV1.vertices;
%  Points2=FV2.vertices;
% % Make binary images form the surface-descriptions
%  V1 = polygon2voxel(FV1,[200 200 150],'none',false);  V1=imfill(V1,'holes');
%  V2 = polygon2voxel(FV2,[200 200 150],'none',false);  V2=imfill(V2,'holes');
% 
% % Do the point to plane matching
%  [Points1n,Points2n]=DistanceFieldCP(Points1,V1,Points2,V2);
%
% % Show the results
%  figure, hold on;
%  plot3(Points1(:,2),Points1(:,1),Points1(:,3),'b.');
%  plot3(Points2(:,2),Points2(:,1),Points2(:,3),'r.');
%  plot3([Points2(:,2) Points2n(:,2)]',[Points2(:,1) Points2n(:,1)]',[Points2(:,3) Points2n(:,3)]','g');
%  plot3([Points1(:,2) Points1n(:,2)]',[Points1(:,1) Points1n(:,1)]',[Points1(:,3) Points1n(:,3)]','m');
%
%
%  Function is written by D.Kroon University of Twente (March 2011)

Points1n=fastICP(I2,I1,Points1);
if(nargout>1)
    Points2n=fastICP(I1,I2,Points2);
end

function Points1n=fastICP(I2,I1,Points1)
s=size(I1);
minp=floor(min(Points1,[],1)-4);
maxp=ceil(max(Points1,[],1)+4);

PadOn=any(minp<1)||any(maxp>s);
if(PadOn)
    % Points close to boundary Padarray
    I1= padarray(I1,[5 5 5],0,'both');
    I2= padarray(I2,[5 5 5],0,'both');
    Points1=Points1+5;
    minp=floor(min(Points1,[],1)-4);
    maxp=ceil(max(Points1,[],1)+4);
    maxp(maxp>s)=s(maxp>s);
    minp(minp<1)=1;
end

I1=I1(minp(1):maxp(1),minp(2):maxp(2),minp(3):maxp(3));
I2=I2(minp(1):maxp(1),minp(2):maxp(2),minp(3):maxp(3));
Points1=bsxfun(@minus,Points1,minp-1);
Points1n=icp(I2,I1,Points1);
Points1n=bsxfun(@plus,Points1n,minp-1);

if(PadOn)
    Points1n=Points1n-5;
end

function Points2n=icp(I1,I2,Points2)
if(size(Points2,2)==2)
    % Calculate the normals of both surfaces
    [ny1,nx1]=gradient(imgaussian(double(I1),1.5));
    [ny2,nx2]=gradient(imgaussian(double(I2),1.5));
    l=sqrt(nx1.^2+ny1.^2)+eps; nx1=nx1./l; ny1=ny1./l;
    l=sqrt(nx2.^2+ny2.^2)+eps; nx2=nx2./l; ny2=ny2./l;
    
    % Get the Normals of the Points form the volume-normals
    nx = interp2(nx2,Points2(:,2),Points2(:,1),'linear')';
    ny = interp2(ny2,Points2(:,2),Points2(:,1),'linear')';
    
    % Get only the voxelized surface of the volume
    B1=xor(I1,imerode(I1,ones(3,3)));
    
    % Divide the normal space in 8 pieces, to match points only with
    % surfaces with approximately the same surface normal
    a=linspace(0,2*pi,8);
    V=[sin(a(:)) cos(a(:))];
    md=V(1,1)*V(2,1)+V(1,2)*V(2,2);
    
    % The closest (interpolated) points
    Points2n=zeros(size(Points2));
    
    % Do the point surface matching for all 8 normal directions
    for i=1:size(V,1)
        % Use the dot product to get only the point and surface pixels
        % with approximately the same normal
        checkb=(nx1.*V(i,1)+ny1*V(i,2))>(md*0.5);
        B1c=B1; B1c(~(checkb))=false;
        check=(nx.*V(i,1)+ny*V(i,2))>=md;
        % Dot the point to plane matching
        Points2n(check,:)=nearestplanepart(B1c,Points2(check,:));
    end
else
    % Calculate the normals of both surfaces
    [ny1,nx1,nz1]=gradient(imgaussian(double(I1),1));
    [ny2,nx2,nz2]=gradient(imgaussian(double(I2),1));
    l=sqrt(nx1.^2+ny1.^2+nz1.^2)+eps; nx1=nx1./l; ny1=ny1./l; nz1=nz1./l;
    l=sqrt(nx2.^2+ny2.^2+nz2.^2)+eps; nx2=nx2./l; ny2=ny2./l; nz2=nz2./l;
    
    % Get the Normals of the Points form the volume-normals
    nx = interp3(nx2,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    ny = interp3(ny2,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    nz = interp3(nz2,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    
    % Detect Points outside the volume
    outside=isnan(nx(:));

    % Get only the voxelized surface of the volume
    B1=xor(I1,imerode(I1,ones(3,3,3)));
    
    % Divide the normal space in 8 pieces, to match points only with
    % surfaces with approximately the same surface normal
    V=[ 0.8507, 0.5257,0; -0.8507, 0.5257,0; -0.8507, -0.5257,0; 0.8507, -0.5257,0;
        0.5257,0, 0.8507; 0.5257,0, -0.8507; -0.5257,0, -0.8507; -0.5257,0, 0.8507;
        0, 0.8507, 0.5257; 0, -0.8507, 0.5257; 0, -0.8507, -0.5257; 0, 0.8507, -0.5257];
    md=V(1,1)*V(2,1) + V(1,2)*V(2,2) + V(1,3)*V(2,3);
    
    % The closest (interpolated) points
    Points2n=zeros(size(Points2));
    
    % Do the point surface matching for all 8 normal directions
    for i=1:size(V,1)
        % Use the dot product to get only the point and surface pixels
        % with approximately the same normal
        checkb=(nx1.*V(i,1)+ny1*V(i,2)+nz1*V(i,3))>(md*0.5);
        B1c=B1; B1c(~(checkb))=false;
        check=(nx.*V(i,1)+ny*V(i,2)+nz*V(i,3))>=md;
        % Dot the point to plane matching
        Points2n(check,:)=nearestplanepart(B1c,Points2(check,:));
    end
    
    % Points outside the volume must stay on the original position
    Points2n(outside,:)=Points2(outside,:);
end

function Points3=nearestplanepart(B1,Points2)
% Make a distance map from the (boundary) object in the logical volume
D1 = double(bwdist(B1,'euclidean'));
if(size(Points2,2)==2)
    % Calculate the closest points from the distance map
    [fy,fx]=gradient(D1);
    d = -interp2(D1,Points2(:,2),Points2(:,1),'linear')';
    dx = interp2(fx,Points2(:,2),Points2(:,1),'linear')';
    dy = interp2(fy,Points2(:,2),Points2(:,1),'linear')';
    Points3=[dx(:).*d(:) dy(:).*d(:)]+Points2;
else
    % Calculate the closest points from the distance map
    [fy,fx,fz]=gradient(D1);
    d = -interp3(D1,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    dx = interp3(fx,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    dy = interp3(fy,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    dz = interp3(fz,Points2(:,2),Points2(:,1),Points2(:,3),'linear')';
    Points3=[dx(:).*d(:) dy(:).*d(:) dz(:).*d(:)]+Points2;
end

