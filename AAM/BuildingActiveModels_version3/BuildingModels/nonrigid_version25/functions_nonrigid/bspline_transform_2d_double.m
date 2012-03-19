function [Iout,Tx,Ty]=bspline_transform_2d_double(Ox,Oy,Iin,dx,dy,mode,ImageSize)
% Bspline transformation grid function
% 
% [Iout,Tx,Ty]=bspline_transform_2d_double(Ox,Oy,Iin,dx,dy,mode,ImageSize)
%
% Inputs,
%   Ox, Oy : are the grid points coordinates
%   Iin : is input image, Iout the transformed output image
%   dx and dy :  are the spacing of the b-spline knots
%   mode: If 0: backward linear interpolation and outside pixels set to nearest pixel
%            1: backward linear interpolation and outside pixels set to zero
%            2: backward cubic interpolation and outsite pixels set to nearest pixel
%            3: backward cubic interpolation and outside pixels set to zero
%            4: forward splatting kernel based interpolation
%   (optional)
%   ImageSize : Image output size
%
% Outputs,
%   Iout: The transformed image
%   Tx: The transformation field in x direction
%   Ty: The transformation field in y direction
%
% This function is an implementation of the b-spline registration
% algorithm in "D. Rueckert et al. : Nonrigid Registration Using Free-Form 
% Deformations: Application to Breast MR Images".
% 
% We used "Fumihiko Ino et al. : a data distrubted parallel algortihm for 
%  nonrigid image registration" for the correct formula's, because 
% (most) other papers contain errors. 
%
% Function is written by D.Kroon University of Twente (June 2009)

if(~exist('ImageSize','var')), ImageSize=size(I); end

% Make all x,y indices
if(mode~=4)
	[x,y]=ndgrid(0:ImageSize(1)-1,0:ImageSize(2)-1);
else
	[x,y]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);
end

% Calulate the transformation of all image coordinates by the b-spline grid
O_trans(:,:,1)=Ox; O_trans(:,:,2)=Oy;
Tlocal=bspline_trans_points_double(O_trans,[dx dy],[x(:) y(:)],false);

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

if(nargout>1)
    % Store transformation fields
    Tx=reshape(Tlocal(:,1),[size(Iin,1) size(Iin,2)])-x;  
    Ty=reshape(Tlocal(:,2),[size(Iin,1) size(Iin,2)])-y; 
end
if(mode==4)
    Iout=image_interpolation_forward(Iin,Tlocal(:,1),Tlocal(:,2),ImageSize);
else
    Iout=image_interpolation_backward(Iin,Tlocal(:,1),Tlocal(:,2),Interpolation,Boundary,ImageSize);
end
