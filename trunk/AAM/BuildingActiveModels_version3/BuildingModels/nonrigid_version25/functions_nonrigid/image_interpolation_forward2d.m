function Iout = image_interpolation_forward2d(Iin,Tlocalx,Tlocaly,ImageSize)
% This function is used to transform an 2D image, in a forwards way with an
% transformation image.
%
%   Iout = image_interpolation_forward2d(Iin,Tlocalx,Tlocaly,ImageSize)
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


% Reshape the transformation (coordinates) to transformation images
Tlocalx=reshape(Tlocalx,[size(Iin,1) size(Iin,2)]);
Tlocaly=reshape(Tlocaly,[size(Iin,1) size(Iin,2)]);
if(nargin<4), ImageSize=[size(Iin,1) size(Iin,2)]; end

Iout=zeros([ImageSize size(Iin,3)]);
Icount=zeros(ImageSize);
for jx=1:size(Tlocalx,1);
    for jy=1:size(Tlocalx,2);
        x2=Tlocalx(jx,jy);
        y2=Tlocaly(jx,jy);
        
        x2r=round(x2);
        y2r=round(y2);
        
        jxmin=max(jx-1,1);
        jymin=max(jy-1,1);
        jxmax=min(jx+1,size(Iin,1));
        jymax=min(jy+1,size(Iin,2));
        
        b1=abs(x2-Tlocalx(jxmin,jy));
        b2=abs(x2-Tlocalx(jxmax,jy));
        b3=abs(y2-Tlocaly(jx,jymin));
        b4=abs(y2-Tlocaly(jx,jymax));
        b5=abs(x2-Tlocalx(jx,jymin));
        b6=abs(x2-Tlocalx(jx,jymax));
        b7=abs(y2-Tlocaly(jxmin,jy));
        b8=abs(y2-Tlocaly(jxmax,jy));
        
        b=ceil(max(1.412*max([b1 b2 b3 b4 b5 b6 b7 b8]),1));
        bc=-(1.0/(b*b))/0.04;
                
        outside=(x2r<=-b)||(x2r>ImageSize(1)+b)||(y2r<=-b)||(y2r>ImageSize(2)+b);
        if(outside), continue, end;
        
        inside=(x2r>b)&&(x2r<ImageSize(1)-b)&&(y2r>b)&&(y2r<ImageSize(2)-b);
        if(inside)
            for ix=-b:b
                for iy=-b:b
                    x3=x2r+ix;
                    y3=y2r+iy;
                    d=bc*((x3-x2)^2+(y3-y2)^2);
                    w= exp(d);
                    Iout(x3,y3,:)=Iout(x3,y3,:)+w*Iin(jx,jy,:);
                    Icount(x3,y3)=Icount(x3,y3)+w;
                end
            end
        else
            for ix=-b:b
                for iy=-b:b
                    x3=x2r+ix;
                    y3=y2r+iy;
                    if((x3>0&&x3<=ImageSize(1))&&(y3>0&&y3<=ImageSize(2)))
                        d=bc*((x3-x2)^2+(y3-y2)^2);
                        w= exp(d);
                        Iout(x3,y3,:)=Iout(x3,y3,:)+w*Iin(jx,jy,:);
                        Icount(x3,y3)=Icount(x3,y3)+w;
                    end
                end
            end
        end
    end
end
for i=1:size(Iin,3)
    Iout(:,:,i)=Iout(:,:,i)./max(Icount,eps(0));
end

