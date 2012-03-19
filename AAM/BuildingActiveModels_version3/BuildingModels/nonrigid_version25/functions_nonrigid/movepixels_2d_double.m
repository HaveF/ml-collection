function Iout=movepixels_2d_double(Iin,Tx,Ty,mode)
% This function movepixels, will translate the pixels of an image
%  according to x and y translation images (bilinear interpolated). 
% 
%  Iout = movepixels_2d_double(I,Tx,Ty,mode,ImageSize);
%
% Inputs;
%   Tx, Ty: The transformation images, describing the
%             (backwards) translation of every pixel in x and y direction.
%   mode: If 0: linear interpolation and outside pixels set to nearest pixel
%            1: linear interpolation and outside pixels set to zero
%            (cubic interpolation only supported by compiled mex file)
%            2: cubic interpolation and outsite pixels set to nearest pixel
%            3: cubic interpolation and outside pixels set to zero
%   (optional)
%   ImageSize : Image output size
%
% Outputs,
%   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (February 2009)
  
% Make all x,y indices
[x,y]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);

% Calculate the Transformed coordinates
Tlocalx = x+Tx;
Tlocaly = y+Ty;

switch(mode)
	case 0
		Interpolation='bilinear';
		Boundary='replicate';
	case 1
		Interpolation='bilinear';
		Boundary='zero';
	case 2
		Interpolation='bicubic';
		Boundary='replicate';
    case 4
	otherwise
		Interpolation='bicubic';
		Boundary='zero';
end

if(mode==4)
    Iout=image_interpolation_forward(Iin,Tlocalx,Tlocaly,ImageSize);
else
    Iout=image_interpolation_backward(Iin,Tlocalx,Tlocaly,Interpolation,Boundary,ImageSize);
end
