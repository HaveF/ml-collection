function [I1,Lines,Vertices,Linesl,Verticesl]=GetImageContour(filename)
% Read the Image
I1=imread(filename);
% Convert it to a binary image
I1=mean(double(I1)/255,3)>0.5;
% Keep only the largest object (remove noise pixels)
[L,N] = bwlabel(I1,4);
s=zeros(1,N); for i=1:N, s(i)=sum(L(:)==i); end; [t,ind]=max(s);
I1=double(L==ind);
% Smooth to get a smoother (less discrete) contour
I1=imgaussian(I1,1.5);
% Get the contour of image
scale=0.1;
[Lines,Vertices]=isocontour(imresize(I1,scale),0.5); Vertices=(Vertices-0.5)*(1/scale)+0.5;
[Linesl,Verticesl]=isocontour(I1,0.5);


