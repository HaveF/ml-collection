function [T,Xreg,X1,X2]=point_registration_diff(sizeI,Xmoving,Xstatic,Options,X1,X2)
% This function creates a 2D or 3D b-spline grid, which transform space to
% fit a set of points Xmoving to set of corresponding points in  Xstatic. 
% Usefull for:
% - For image-registration based on corresponding landmarks like in
% Sift or OpenSurf (see Mathworks). 
% - 2D and 3D Spline based Data gridding and surface fitting
% - Smooth filtering of 2d / 3D Point data.
%
%   [T,Xreg]=point_registration(sizeI,Xstatic,Xmoving,Options);
%
% Inputs,
%   sizeI : The size of the (virtual) image/space which will be warped
%            With the b-spline grid
%   Xmoving : List with 2D or 3D points N x 2, or N x 3, these points will be
%            warped be the fitted b-sline grid to transform to the 
%            static points in Xstatic
%   Xstatic : List with 2D or 3D points N x 2, or N x 3, corresponding with Xmoving
%   Options : Struct with options, see below.
%
% Outputs,
%   T : The diffeomorphic transformation fields
%   Xreg: The points in Xmoving transformed by the fitted b-spline grid.
%
% Options,
%   Options.Verbose: Display Debug information 0,1 or 2
%   Options.MaxRef : Maximum number of grid refinements steps (default 5)
%   Options.Forwards : Use (slower) forward interpolation (default false).
%
%
 % Load corresponding landmarks
%      load('images/starpoints.mat');
%    % Load the images
%      I1=im2double(imread('images/star1.png')); 
%      I2=im2double(imread('images/star2.png')); 
%    % Fit the bspline grid to the corresponding landmarks
%      options.Verbose=true;
%      T=point_registration_diff(size(I1),[x1(:) y1(:)],[x2(:) y2(:)],options);
%    % Transform the 2D image  
%      Ireg=movepixels(I1,T,4);
%    % Show the result
%      figure,
%      subplot(2,2,1),imshow(I1); title('Moving Image'); 
%      hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end
%      subplot(2,2,2),imshow(I2); title('Static Image');
%      subplot(2,2,3),imshow(Ireg); title('Registered Image')
%      hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end 
%      Ireg2=movepixels(I2,T,3);
%      subplot(2,2,4),imshow(Ireg2); title('Registered Image')
%      hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end 
% 
%    % Show b-spline grid
%      Igrid=make_grid_image([8 8],size(I1));
%      figure, 
%      subplot(1,2,1), imshow(Igrid)
%      Ireg=movepixels(Igrid,T,3);
%      subplot(1,2,2), imshow(Ireg)
%      hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end
%      
%      
%      
%     I1=im2double(imread('d:\matlab\lena.jpg')); 
%        
%        Xstatic=[1 1;
%        1 256;
%        1 128;
%        256 128;       
%        128-32 128;
%        128+32 128;
%        256 1;
%        256 256];
%   
%     Xmoving=[1 1;
%        1 256;
%        1 128;
%        256 128;
%        128+62 128;
%        128-62 128;
%        256 1;
%        256 256];
%    
%         
%      T=point_registration_diff(size(I1),Xstatic,Xmoving,options);
%       Ireg=movepixels(I1,T,4);
%       figure, 
%       subplot(1,2,1), imshow(I1);
%       subplot(1,2,2), imshow(Ireg);
%
%
%   Xstatic=[1 1 64;
%        1 128 64;
%        64+32 64 64
%        64-32 64 64
%        128 1 64;
%        128 128 64];
%   
%     Xmoving=[1 1 64;
%        1 128 64;
%        64-32 64 64
%        64+32 64 64
%        128 1 64;
%        128 128 64];
%     option=struct; options.MaxRef=3;
%      options.Verbose=true;
%     sizeI=[128 128 128]; 
%     T=point_registration_diff(sizeI,Xstatic,Xmoving,options);
%  Igrid=make_grid_image([8 8 8],sizeI);
%  Ireg=movepixels(Igrid,T,4);
% showcs3(Ireg);
    
%  Function is written by D.Kroon University of Twente (May 2011)
  
% add all needed function paths
add_function_paths;

% Disable warning
warning('off', 'MATLAB:maxNumCompThreads:Deprecated')

% Process inputs
defaultoptions=struct('Verbose',false,'MaxRef',5,'Forwards',false);
if(~exist('Options','var')), Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags), if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end, end
    if(length(tags)~=length(fieldnames(Options))),
        warning('register_images:unknownoption','unknown options found');
    end
end

% Check if all points are inside the b-spline volume
outside=false;
if(any(Xmoving(:)<1)||any(Xmoving(:,1)>sizeI(1))||any(Xmoving(:,2)>sizeI(2))), outside=true; end
if(any(Xstatic(:)<1)||any(Xstatic(:,1)>sizeI(1))||any(Xstatic(:,2)>sizeI(2))), outside=true; end
switch(size(Xstatic,2))
    case 3
		% Check if all points are inside the b-spline volume
		if(any(Xmoving(:,3)>sizeI(3))), outside=true; end
		if(any(Xstatic(:,3)>sizeI(3))), outside=true; end
end        

if(outside)
    warning('point_registration_diff:boundary','There are points outside the image boundary defined by Isize');
    sizeI=ceil(max([sizeI;max(Xmoving,[],1);max(Xstatic,[],1)],[],1));
end

switch(size(Xstatic,2))
    case 2
        % Calculate max refinements steps
        MaxItt=min(floor(log2(sizeI(1:2)/2)));

        % set b-spline grid spacing in x and y direction
        Spacing=[2^MaxItt 2^MaxItt];
    case 3
	       % Calculate max refinements steps
        MaxItt=min(floor(log2(sizeI(1:3)/2)));

        % set b-spline grid spacing in x,y and z direction
        Spacing=[2^MaxItt 2^MaxItt 2^MaxItt];
end


% Make an initial uniform b-spline grid
O_ref = make_init_grid(Spacing,sizeI);

% Calculate difference between the points
R=Xstatic-Xmoving;

switch(size(Xstatic,2))
    case 2
        T=zeros([sizeI(1) sizeI(2) 2]);
        [xg,yg]=ndgrid(1:sizeI(1),1:sizeI(2));
        T(:,:,1)=xg;
        T(:,:,2)=yg;
    case 3
        T=zeros([sizeI(1) sizeI(2) sizeI(3) 3],'single');
        [xg,yg,zg]=ndgrid(1:sizeI(1),1:sizeI(2),1:sizeI(3));
        T(:,:,:,1)=xg;
        T(:,:,:,2)=yg;
        T(:,:,:,3)=zg;
end

Xreg=Xmoving;

% Loop through all refinement itterations
%err=1e200;
for i=1:Options.MaxRef
    for k=1:8
        if(Options.Verbose)
            disp('.');
            disp(['Iteration : ' num2str(i) '/' num2str(Options.MaxRef)]);
            disp(['Grid size : ',num2str(size(O_ref))]);
        end

        % Make a b-spline grid which minimizes the difference between the
        % corresponding points
        O=bspline_grid_fitting(zeros(size(O_ref)),Spacing,R,Xreg);
        c=min(Spacing)*0.35;

        O(O>c )= c;
        O(O<-c)=-c;
        O=O+O_ref;
        
        switch(size(Xstatic,2))
            case 2
                Tx=T(:,:,1); 
                Ty=T(:,:,2);
                Tuxy=bspline_trans_points_double(O,Spacing,[Tx(:) Ty(:)]);
                Xreg=bspline_trans_points_double(O,Spacing,Xreg);
                if(nargin>4)
                    X1=bspline_trans_points_double(O,Spacing,X1);
                end
                if(nargin>5)
                    X2=bspline_trans_points_double(O,Spacing,X2);
                end
                Tu(:,:,1)=reshape(Tuxy(:,1),[sizeI(1) sizeI(2)]);
                Tu(:,:,2)=reshape(Tuxy(:,2),[sizeI(1) sizeI(2)]);
            case 3
                Tx=T(:,:,:,1); 
                Ty=T(:,:,:,2);
                Tz=T(:,:,:,3);
                Tuxy=bspline_trans_points_double(O,Spacing,[Tx(:) Ty(:) Tz(:)]);
                Xreg=bspline_trans_points_double(O,Spacing,Xreg);
                if(nargin>4)
                    X1=bspline_trans_points_double(O,Spacing,X1);
                end
                if(nargin>5)    
                    X2=bspline_trans_points_double(O,Spacing,X2);
                end
                Tu(:,:,:,1)=reshape(Tuxy(:,1),[sizeI(1) sizeI(2) sizeI(3)]);
                Tu(:,:,:,2)=reshape(Tuxy(:,2),[sizeI(1) sizeI(2) sizeI(3)]);
                Tu(:,:,:,3)=reshape(Tuxy(:,3),[sizeI(1) sizeI(2) sizeI(3)]);
        end
        T=Tu;

        % Calculate the remaining difference between the points
        R=Xstatic-Xreg;
        %errold=err; 
        err=mean(sqrt(sum(R.^2,2)));
        if(Options.Verbose)
            disp(['Mean Distance : ',num2str(err)]);
        end
        %if((err/errold)>0.98), break; end
    end
    
    if(i<Options.MaxRef)
        % Refine the update-grid and reference grid
        switch(size(Xstatic,2))
            case 2
                [O_ref ,Spacing]=refine_grid(O_ref,Spacing,sizeI(1:2));
            case 3
                [O_ref ,Spacing]=refine_grid(O_ref,Spacing,sizeI);
        end        
    end
end
switch(size(Xstatic,2))
    case 2
        T(:,:,1)= T(:,:,1)-xg;
        T(:,:,2)= T(:,:,2)-yg;
    case 3
        T(:,:,:,1)=T(:,:,:,1)-xg;
        T(:,:,:,2)=T(:,:,:,2)-yg;
        T(:,:,:,3)=T(:,:,:,3)-zg;
end

% The final transformation grid, is the reference grid with update-grid
% added

% Check if Jacobian > EPS (no folding or mesh collapse)
Err=DiscreteJacobian(T);
if(Err)
    switch(size(Xstatic,2))
    case 2
        b=ceil(min([size(T,1) size(T,2)])/10);
        Err=DiscreteJacobian(T(b+1:end-b,b+1:end-b,:));
    case 3
        b=ceil(min([size(T,1) size(T,2) size(T,3)])/10);
        Err=DiscreteJacobian(T(b+1:end-b,b+1:end-b,b+1:end-b,:));
    end        
    if(Err)
        warning('point_registration_diff:jacobian','Jacobian zero or negative in center of image');
    else
        disp('Jacobian zero or negative at boundary');
    end
end


function add_function_paths()
try
    functionname='point_registration.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions'])
    addpath([functiondir '/functions_affine'])
    addpath([functiondir '/functions_nonrigid'])
catch me
    disp(me.message);
end

 