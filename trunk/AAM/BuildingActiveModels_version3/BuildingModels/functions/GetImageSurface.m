function [I1,Faces,Vertices,Facesl,Verticesl]=GetImageSurface(filename)
% Read the Image - argument list modified
d=load(filename);

% Keep only the largest object (remove noise pixels)
[L,N] = bwlabeln(d.V,6);
s=zeros(1,N); for i=1:N, s(i)=sum(L(:)==i); end; [t,ind]=max(s);
I1=single(L==ind);

% Smooth to get a smoother (less discrete) contour

I1s=imgaussian(I1,1.5);
% Get the contour of image-volume

scale=1/4;
Surface_JAW_Small=isosurface(padarray(imresize3d(single(I1s),scale,[],'linear'),[1 1 1],0,'both'),0.4);
Surface_JAW_Small.vertices=Surface_JAW_Small.vertices-1;
Surface_JAW_Small.vertices=(1/scale)*(Surface_JAW_Small.vertices(:,[2 1 3]));
Vertices=Surface_JAW_Small.vertices;
Faces=Surface_JAW_Small.faces;

scale=1/2;
Surface_JAW_Large=isosurface(padarray(imresize3d(single(I1s),scale,[],'linear'),[1 1 1],0,'both'),0.4);
Surface_JAW_Large.vertices=Surface_JAW_Large.vertices-1;
Surface_JAW_Large.vertices=(1/scale)*(Surface_JAW_Large.vertices(:,[2 1 3]));
Facesl=Surface_JAW_Large.faces;
Verticesl=Surface_JAW_Large.vertices;
