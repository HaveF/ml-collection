function q_u=flbp81_rgb_PDF(Image,center,w_halfsize,lbpThreshold,redBins,greenBins,blueBins)
% *************************************************************************
% function: calculate the target model with color-texture histogram
% parameters: 
%       Image: the 1th image
%       center: the center of target window
%       w_halfsize: the bandwidth of target window
%       lbpThreshold,
%       redBins,greenBins,blueBins quantification scheme of rgb space
% returns:
%       q_u: target model
%
% Authors: Jifeng Ning£¬Lei Zhang, David Zhang, Chengke Wu
% *************************************************************************

grayImg=rgb2gray(uint8(Image));                     % convert rgb image into gray image
[height,width]=size(grayImg);
tempImg=zeros(height+2,width+2);
tempImg(2:height+1,2:width+1)=grayImg;              % 

textureMat=LBPTexture(tempImg,8,1,lbpThreshold); % compute lbp texture image with lbp81 operator
% preserve the texture pattern related to line and corner point
textureMat(textureMat==1 | textureMat==2 | textureMat==8 | textureMat==9 | textureMat==10)=0; 
textureMat(textureMat>0)=textureMat(textureMat>0)-2;

sum_q=0;
histo=zeros(redBins,greenBins,blueBins,5);         % initialize the target histogram

rmin=center(1)-w_halfsize(1);
rmax=center(1)+w_halfsize(1);
cmin=center(2)-w_halfsize(2);
cmax=center(2)+w_halfsize(2);


wmax=(rmin-center(1)).^2+(cmin-center(2)).^2+1;
for i=rmin:rmax  % 
    for j=cmin:cmax
        if textureMat(i,j)~=0
            d=(i-center(1)).^2+(j-center(2)).^2;
            w=wmax-d;  % weight value corresponding to pixel(i,j)
            R=floor(Image(i,j,1)/(256/redBins))+1;
            G=floor(Image(i,j,2)/(256/greenBins))+1;
            B=floor(Image(i,j,3)/(256/blueBins))+1;
            T=textureMat(i,j);
            histo(R,G,B,T)=histo(R,G,B,T)+w;  % color-texture joint histogram
        end
    end
end

% convert histo (matrix) into 1D vector
for i=1:redBins
    for j=1:greenBins
        for k=1:blueBins
            for s=1:5
                index=(i-1)*greenBins*blueBins*5+(j-1)*blueBins*5+(k-1)*5+s;
                q_u(index)=histo(i,j,k,s);
                sum_q=sum_q+q_u(index);
            end
        end
    end
end

q_u=q_u/sum_q; 
return;