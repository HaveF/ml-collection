function Iout = image_interpolation_forward(Iin,Tlocalx,Tlocaly,ImageSize)
% This function is used to transform an 2D image, in a forwards way with an
% transformation image.
%
%   Iout = image_interpolation_forward(Iin,Tlocalx,Tlocaly,ImageSize)
%
% inputs,
%	   Iin : 2D greyscale or color input image
%	   Tlocalx,Tlocaly : (Forwards) Transformation images for all image pixels
%
%	(optional)
%	   ImageSize:    - Size of output image
% outputs,
%  	   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (September 2010)
%
% % Example
%   I=im2double(imread('d:\matlab\lena.jpg'));
%   [x,y]=ndgrid(1:size(I,1),1:size(I,2));
%   ImageSize=[256 256];
%   Tlocalx=(x-128)*2-(y-128)*2;
%   Tlocaly=(y-128)*2+(x-128)*2;
%   J=image_interpolation_forward(I,Tlocalx,Tlocaly,ImageSize);
%   figure, imshow(J);
%
% % Example 3D
%   I=make_grid_image([16,16,16],[256 256 256]);
%   [x,y,z]=ndgrid(1:size(I,1),1:size(I,2),1:size(I,3));
%   ImageSize=[256 256 256];
%   Tlocalx=x; %(x-128)*2-(y-128)*2;
%   Tlocaly=y; %(y-128)*2+(x-128)*2;
%   Tlocalz=z; %(y-128)*2+(x-128)*2;
%   J=image_interpolation_forward3d(I,Tlocalx,Tlocaly,Tlocalz,ImageSize);

% Reshape the transformation (coordinates) to transformation images
Tlocalx=reshape(Tlocalx,[size(Iin,1) size(Iin,2)]);
Tlocaly=reshape(Tlocaly,[size(Iin,1) size(Iin,2)]);
if(nargin<4), ImageSize=[size(Iin,1) size(Iin,2)]; end

Iout = image_interpolation_forward2d(Iin,Tlocalx,Tlocaly,ImageSize);