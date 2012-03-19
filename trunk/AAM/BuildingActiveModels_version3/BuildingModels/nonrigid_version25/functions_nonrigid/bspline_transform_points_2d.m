function Tlocal=bspline_transform_points_2d(O_trans,Spacing,X)

% Make row vectors of input coordinates
x2=X(:,1); y2=X(:,2);

% This code calculates for every coordinate in X, the indices of all
% b-spline knots which have influence on the transformation value of
% this point
[m,l]=ndgrid(0:3,0:3); m=m(:)'; l=l(:)';
ixs=floor(x2/Spacing(1));
iys=floor(y2/Spacing(2));
ix=repmat(ixs,[1 16])+repmat(m,[length(x2) 1]); ix=ix(:);
iy=repmat(iys,[1 16])+repmat(l,[length(y2) 1]); iy=iy(:);

% Size of the b-spline grid
s=size(O_trans);

% Points outside the bspline grid are set to the upper corner
Check_bound=(ix<0)|(ix>(s(1)-1))|(iy<0)|(iy>(s(2)-1));
ix(Check_bound)=1; iy(Check_bound)=1;
Check_bound_inv=double(~Check_bound);

% Look up the b-spline knot values in neighborhood of the points in (x2,y2)
Cx=O_trans(ix+iy*s(1)+1).*Check_bound_inv;
Cx=reshape(Cx,[length(x2) 16]);
Cy=O_trans(ix+iy*s(1)+s(1)*s(2)+1).*Check_bound_inv;
Cy=reshape(Cy,[length(x2) 16]);

% Calculate the b-spline interpolation constants u,v in the center cell
% range between 0 and 1
v  = (x2-ixs*Spacing(1))/Spacing(1);
u  = (y2-iys*Spacing(2))/Spacing(2);

% Get the b-spline coefficients in a matrix W, which contains
% the influence of all knots on the points in (x2,y2)
W=bspline_coefficients(v,u);

% Calculate the transformation of the points in (x2,y2) by the b-spline grid
Tlocal(:,1)=sum(W.*Cx,2);
Tlocal(:,2)=sum(W.*Cy,2);

