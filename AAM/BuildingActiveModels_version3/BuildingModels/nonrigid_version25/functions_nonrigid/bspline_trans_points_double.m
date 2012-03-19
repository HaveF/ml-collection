function [Tlocal,Dlocal]=bspline_trans_points_double(O_trans,Spacing,X,check)
% If Spacing and points are not an integer, we can not use fast look up tables,
% but have to calculate the bspline coefficients for every point
if(nargin<4)
    check = any(mod(Spacing,1)>0)|any(mod(X(:),1)>0);
end
O_trans=double(O_trans); Spacing=double(Spacing);
switch size(X,2)
    case 2,
        if(check)
            Tlocal=bspline_transform_points_2d(O_trans,Spacing,X);
        else
            if(nargout>1)
                [Tlocal,Dlocal]=bspline_transform_fast_2d(O_trans,Spacing,X);
            else
                Tlocal=bspline_transform_fast_2d(O_trans,Spacing,X);
            end
        end
    case 3,
    Tlocal=zeros(size(X),class(X));
    ind=round(linspace(1,size(X,1),10));
    for q=1:(length(ind)-1),
        Tlocal(ind(q):ind(q+1),:)=bspline_transform_points_3d(O_trans,Spacing,double(X(ind(q):ind(q+1),:)));
    end
end

function [Tlocal,Dlocal]=bspline_transform_fast_2d(O_trans,Spacing,X)
% Make row vectors of input coordinates
x2=X(:,1); y2=X(:,2);

% Make polynomial look up tables
Bu=zeros(4,Spacing(1));
Bv=zeros(4,Spacing(2));
Bdu=zeros(4,Spacing(1));
Bdv=zeros(4,Spacing(2));

x=0:Spacing(1)-1;
u=(x/Spacing(1))-floor(x/Spacing(1));
Bu(0*Spacing(1)+x+1) = (1-u).^3/6;
Bu(1*Spacing(1)+x+1) = ( 3*u.^3 - 6*u.^2 + 4)/6;
Bu(2*Spacing(1)+x+1) = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
Bu(3*Spacing(1)+x+1) = u.^3/6;

y=0:Spacing(2)-1;
v=(y/Spacing(2))-floor(y/Spacing(2));
Bv(0*Spacing(2)+y+1) = (1-v).^3/6;
Bv(1*Spacing(2)+y+1) = ( 3*v.^3 - 6*v.^2 + 4)/6;
Bv(2*Spacing(2)+y+1) = (-3*v.^3 + 3*v.^2 + 3*v + 1)/6;
Bv(3*Spacing(2)+y+1) = v.^3/6;

if(nargout>1)
    Bdu(0*Spacing(1)+x+1) = -(u - 1).^2/2;
    Bdu(1*Spacing(1)+x+1) = (3*u.^2)/2 - 2*u;
    Bdu(2*Spacing(1)+x+1) = -(3*u.^2)/2 + u + 1/2;
    Bdu(3*Spacing(1)+x+1) = u.^2/2;
    
    Bdv(0*Spacing(2)+y+1) = -(v - 1).^2/2;
    Bdv(1*Spacing(2)+y+1) = (3*v.^2)/2 - 2*v;
    Bdv(2*Spacing(2)+y+1) = -(3*v.^2)/2 + v + 1/2;
    Bdv(3*Spacing(2)+y+1) = v.^2/2;
    Bdu=Bdu./Spacing(1);
    Bdv=Bdv./Spacing(2);
end

% Calculate the indexes need to loop up the B-spline values.
u_index=mod(x2,Spacing(1));
v_index=mod(y2,Spacing(2));

i=floor(x2/Spacing(1)); % (first row outside image against boundary artefacts)
j=floor(y2/Spacing(2));

% This part calculates the coordinates of the pixel
% which will be transformed to the current x,y pixel.

Ox=O_trans(:,:,1);
Oy=O_trans(:,:,2);
Tlocalx=0; Tlocaly=0;
Tlocaldxx=0;
Tlocaldyy=0;
Tlocaldxy=0;
Tlocaldyx=0;

a=zeros(size(X,1),4);
b=zeros(size(X,1),4);
if(nargout>1)
    ad=zeros(size(X,1),4);
    bd=zeros(size(X,1),4);
end

IndexO1l=zeros(size(X,1),4);
IndexO2l=zeros(size(X,1),4);
Check_bound1=false(size(X,1),4);
Check_bound2=false(size(X,1),4);

for r=0:3,
    a(:,r+1)=Bu(r*Spacing(1)+u_index(:)+1);
    b(:,r+1)=Bv(r*Spacing(2)+v_index(:)+1);
    if(nargout>1)
        ad(:,r+1)=Bdu(r*Spacing(1)+u_index(:)+1);
        bd(:,r+1)=Bdv(r*Spacing(2)+v_index(:)+1);
    end
    IndexO1l(:,r+1)=(i+r);
    IndexO2l(:,r+1)=(j+r);
    Check_bound1(:,r+1)=(IndexO1l(:,r+1)<0)|(IndexO1l(:,r+1)>(size(O_trans,1)-1));
    Check_bound2(:,r+1)=(IndexO2l(:,r+1)<0)|(IndexO2l(:,r+1)>(size(O_trans,2)-1));
end

for l=0:3,
    for m=0:3,
        IndexO1= IndexO1l(:,l+1);
        IndexO2= IndexO2l(:,m+1);
        Check_bound=Check_bound1(:,l+1)|Check_bound1(:,m+1);
        IndexO1(Check_bound)=1;
        IndexO2(Check_bound)=1;
        Check_bound_inv=double(~Check_bound);
        
        ab=a(:,l+1).*b(:,m+1);
        if(nargout>1)
            abx=ad(:,l+1).*b(:,m+1);
            aby=a(:,l+1).*bd(:,m+1);
        end
        
        c=Ox(IndexO1(:)+IndexO2(:)*size(Ox,1)+1);
        Tlocalx=Tlocalx+Check_bound_inv(:).*ab.*c;
        if(nargout>1)
            Tlocaldxx=Tlocaldxx+Check_bound_inv(:).*abx.*c;
            Tlocaldxy=Tlocaldxy+Check_bound_inv(:).*aby.*c;
        end
        
        c=Oy(IndexO1(:)+IndexO2(:)*size(Oy,1)+1);
        Tlocaly=Tlocaly+Check_bound_inv(:).*ab.*c;
        if(nargout>1)
            Tlocaldyx=Tlocaldyx+Check_bound_inv(:).*abx.*c;
            Tlocaldyy=Tlocaldyy+Check_bound_inv(:).*aby.*c;
        end
    end
end
Tlocal(:,1)=Tlocalx(:);
Tlocal(:,2)=Tlocaly(:);
if(nargout>1)
    Dlocal=(Tlocaldxx).*(Tlocaldyy)-Tlocaldyx.*Tlocaldxy;
end

