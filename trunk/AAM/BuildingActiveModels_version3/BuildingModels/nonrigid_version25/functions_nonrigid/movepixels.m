function I3=movepixels(I1,T,mode,ImageSize)
% This function movepixels, will (backwards) translate the pixels 
% of an 2D/3D image according to x, y (and z) translation images 
% (bilinear interpolated).
% The function is a wrapper around mex files movepixels_2d_double.c and
% movepixels_3d_double.c and movepixels_3d_single.c
%
% J = movepixels(I,T,mode,ImageSize);
%
% Inputs;
%   T : The transformation images, describing the
%            translation of every pixel in x,y and z direction.
%   mode: If 0: backward linear interpolation and outside pixels set to nearest pixel
%            1: backward linear interpolation and outside pixels set to zero
%            2: backward cubic interpolation and outsite pixels set to nearest pixel
%            3: backward cubic interpolation and outside pixels set to zero
%            4: forward splatting kernel based interpolation
%   (optional)
%   ImageSize : Image output size
%
% Outputs,
%   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (March 2009)

if(~exist('mode','var')), mode=0; end
if(~exist('ImageSize','var')), ImageSize=size(I1); end

if(size(I1,3)<4)
    I3=movepixels_2d_double(double(I1),double(T(:,:,1)),double(T(:,:,2)),double(mode),double(ImageSize));
else
    if(isa(I1,'double'))
        I3=movepixels_3d_double(double(I1),double(T(:,:,:,1)),double(T(:,:,:,2)),double(T(:,:,:,3)),double(mode),double(ImageSize));
    else
        I3=movepixels_3d_single(single(I1),single(T(:,:,:,1)),single(T(:,:,:,2)),single(T(:,:,:,3)),single(mode),single(ImageSize));
    end
end
if(~isa(I1,'double')&&~isa(I1,'single'))
	I3=cast(I3,class(I1));
end


