function [O_transO,Spacing]=MakeDiffeomorphic(O_trans,Spacing,sizeI)
% This function MakeDiffeomorphic will make the b-spline grid diffeomorphic
% thus will regularize the grid to have a Jacobian larger then zero
% for the whole b-spline grid.
%
% [Grid,Spacing]=point_registration_diffeomorphic(Grid,Spacing,sizeI)
%
%  Grid: The b-spline controlpoints, can be used to transform another
%        image in the same way: I=bspline_transform(Grid,I,Spacing);
%  Spacing: The uniform b-spline knot spacing
%  sizeI : Size of the image (volume)
%
% Example, 2D Diffeomorphic Warp
%    Xstatic=[1 1;
%      1 128;
%      64+32 64
%      64-32 64
%      128 1;
%      128 128];
% 
%   Xmoving=[1 1;
%      1 128;
%      64-32 64
%      64+32 64
%      128 1;
%      128 128];
%   option=struct; options.MaxRef=4;
%   sizeI=[128 128]; 
%   [O_trans,Spacing]=point_registration(sizeI,Xstatic,Xmoving,options);
%   [O_trans,Spacing]=MakeDiffeomorphic(O_trans,Spacing,sizeI);
%
%   Igrid=make_grid_image(Spacing*2,sizeI);
%   [Ireg,B]=bspline_transform(O_trans,Igrid,Spacing,3);
%   figure, imshow(Ireg)
%
% Example, 3D Diffeomorphic Warp
%   Xstatic=[1 1 64;
%      1 128 64;
%      64+32 64 64
%      64-32 64 64
%      128 1 64;
%      128 128 64];
% 
%   Xmoving=[1 1 64;
%      1 128 64;
%      64-32 64 64
%      64+32 64 64
%      128 1 64;
%      128 128 64];
%   option=struct; options.MaxRef=3;
%   sizeI=[128 128 128]; 
%   [O_trans,Spacing]=point_registration(sizeI,Xstatic,Xmoving,options);
%   [O_trans,Spacing]=MakeDiffeomorphic(O_trans,Spacing,sizeI);
%
%   Igrid=make_grid_image(Spacing,sizeI);
%   [Ireg,B]=bspline_transform(O_trans,Igrid,Spacing,3);
%   showcs3(Ireg);
%
%   Function is written by D.Kroon University of Twente (March 2011)

% Make Diffeomorphic
O_transO=O_trans;
optim=struct('GradObj','on','rho',0.01,'sigma',1,'GoalsExactAchieve',1,'StoreN',10,'HessUpdate','lbfgs','Display','off','MaxIter',100,'DiffMinChange',1e-4,'DiffMaxChange',0.1,'MaxFunEvals',1000,'TolX',0.01,'TolFun',1e-8);
c=[10 30 90 270];
for i=1:length(c);
    ne=CheckJacobianErrors(O_transO,Spacing,sizeI);
    if(ne==0), return; end
    O_transO = fminlbfgs(@(x)jacobiandet_cost_gradient(x,size(O_trans),O_trans,sizeI,Spacing,c(i)),O_trans,optim);
end
ne=CheckJacobianErrors(O_transO,Spacing,sizeI);
if(ne==0), return; end
warning('MakeDiffeomorphic:jacobian', ['jacobian error pixels : ' num2str(ne)]);

function ne=CheckJacobianErrors(O_trans,Spacing,sizeI)
if(size(O_trans,4)>1)
      D=jacobiandet_transform_3d_double(O_trans(:,:,:,1),O_trans(:,:,:,2),O_trans(:,:,:,3),sizeI,Spacing(1),Spacing(2),Spacing(3));
else
      D=jacobiandet_transform_2d_double(O_trans(:,:,1),O_trans(:,:,2),sizeI,Spacing(1),Spacing(2));
end
ne=sum((D(:)<eps));
