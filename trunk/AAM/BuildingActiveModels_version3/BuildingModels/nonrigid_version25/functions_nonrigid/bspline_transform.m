function [J,T]=bspline_transform(O,I,Spacing,mode,ImageSize)
% Function bspline_transform, is a wrapper of the mex 
% bspline_transform_2d_double and bspline_transform_3d mex functions
%
% [Iout,T] = bspline_transform(O,Iin,Spacing,mode,ImageSize)
%
% inputs,
%   Iin :  Input image.
%   O  : Transformation grid of control points
%   Spacing : Are the b-spline grid knot spacings.
%   mode: If 0: backward linear interpolation and outside pixels set to nearest pixel
%            1: backward linear interpolation and outside pixels set to zero
%            2: backward cubic interpolation and outsite pixels set to nearest pixel
%            3: backward cubic interpolation and outside pixels set to zero
%            4: forward splatting kernel based interpolation
%   (optional)
%   ImageSize : Image output size
%
% outputs,
%   Iout : The Rueckert transformed image
%   T : The transformation images (in 2D Tx=T(:,:,1)), describing the
%             (backwards) translation of every pixel in x,y and z direction.
%
% Function is written by D.Kroon University of Twente (February 2009)

% Check if spacing has integer values
if(sum(Spacing-floor(Spacing))>0), error('Spacing must be a integer'); end
if(~exist('mode','var')), mode=0; end
if(~exist('ImageSize','var')), ImageSize=size(I); end

if(size(I,3)<4)
    if(nargout > 1 )
        [J,Tx,Ty]=bspline_transform_2d_double(double(O(:,:,1)),double(O(:,:,2)),double(I),double(Spacing(1)),double(Spacing(2)),double(mode),double(ImageSize));
        T(:,:,1)=Tx; 
        T(:,:,2)=Ty;
    else
        J=bspline_transform_2d_double(double(O(:,:,1)),double(O(:,:,2)),double(I),double(Spacing(1)),double(Spacing(2)),double(mode),double(ImageSize));
        T=0;
    end
	if(~isa(I,'double')), J=cast(J,class(I)); end

else
    if(isa(I,'double'))
        if(nargout > 1 )
            [J,Tx,Ty,Tz]=bspline_transform_3d_double(double(O(:,:,:,1)),double(O(:,:,:,2)),double(O(:,:,:,3)),double(I),double(Spacing(1)),double(Spacing(2)),double(Spacing(3)),double(mode),double(ImageSize));
            T(:,:,:,1)=Tx; T(:,:,:,2)=Ty; T(:,:,:,3)=Tz;
        else
            J=bspline_transform_3d_double(double(O(:,:,:,1)),double(O(:,:,:,2)),double(O(:,:,:,3)),double(I),double(Spacing(1)),double(Spacing(2)),double(Spacing(3)),double(mode),double(ImageSize));
            T=0;
        end
    else
        if(nargout > 1 )
             [J,Tx,Ty,Tz]=bspline_transform_3d_single(single(O(:,:,:,1)),single(O(:,:,:,2)),single(O(:,:,:,3)),single(I),single(Spacing(1)),single(Spacing(2)),single(Spacing(3)),single(mode),single(ImageSize));
             T(:,:,:,1)=Tx; T(:,:,:,2)=Ty; T(:,:,:,3)=Tz;
        else
             J=bspline_transform_3d_single(single(O(:,:,:,1)),single(O(:,:,:,2)),single(O(:,:,:,3)),single(I),single(Spacing(1)),single(Spacing(2)),single(Spacing(3)),single(mode),single(ImageSize));
			T=0;
        end
		J=cast(J,class(I));
	end
end



