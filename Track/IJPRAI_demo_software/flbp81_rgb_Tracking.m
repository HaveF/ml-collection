function [center,iterations]=flbp81_rgb_Tracking(Image,center,w_halfsize,q_u,redBins,greenBins,blueBins,minDist,maxIterNum,incre,lbpThreshold)
% ********************************************************************
% function: mean shift tracking with color-texture joint histogram
% parameters:
%       Image         current frame to track
%       center        initial location of target
%       w_halfsizet   bandwidth of target window
%       q_u           target model
%       redBins,greenBins,blueBins: quantification scheme of rgb space
%       minDist       convergence condition
%       maxIterNum    maxmimal iterative number
%       incre         increase the candidate winodw to robust tracking
%       lbpThreshold  
% return:
%       center        center of tracking-result winodw
% ****************************************************************

grayImg=rgb2gray(uint8(Image));                             % convert rgb image into gray image
[height,width]=size(grayImg);
tempImg=zeros(height+2,width+2);
tempImg(2:height+1,2:width+1)=grayImg;                      

textureMat=LBPTexture(tempImg,8,1,lbpThreshold);           % compute lbp texture image with lbp81 operator

% preserve the texture pattern related to line and corner point
textureMat(textureMat==1 | textureMat==2 | textureMat==8 | textureMat==9 | textureMat==10)=0; % 只保留与边界有关的纹理 模式(1和9代表图像中的光滑区域)
textureMat(textureMat>0)=textureMat(textureMat>0)-2;

sum_p=0; 
histo=zeros(redBins,greenBins,blueBins,5);         % initialize the candidate model histogram

iterations=0;                                      % initialize iteration number
center_old=center;

rmin=center(1)-w_halfsize(1)-incre;
rmax=center(1)+w_halfsize(1)+incre;
cmin=center(2)-w_halfsize(2)-incre;
cmax=center(2)+w_halfsize(2)+incre;
  
while 1
    wmax=(rmin-center(1)).^2+(cmin-center(2)).^2+1;
    for i=rmin:rmax                                % calculate the candidate model
        for j=cmin:cmax
            if textureMat(i,j)~=0
                if (i>=1 & i<=height & j>=1 & j<=width)
                    d=(i-center(1)).^2+(j-center(2)).^2;
                    w=wmax-d;           % weight value corresponding to pixel(i,j)
                    R=floor(Image(i,j,1)/(256/redBins))+1;
                    G=floor(Image(i,j,2)/(256/greenBins))+1;
                    B=floor(Image(i,j,3)/(256/blueBins))+1;
                    T=textureMat(i,j);
                    histo(R,G,B,T)=histo(R,G,B,T)+w;  % color-texture joint histogram
                end
            end
        end
    end

    for i=1:redBins
        for j=1:greenBins
            for k=1:blueBins
                for s=1:5
                    index=(i-1)*greenBins*blueBins*5+(j-1)*blueBins*5+(k-1)*5+s;
                    p_u(index)=histo(i,j,k,s);
                    sum_p=sum_p+p_u(index);
                end
            end
        end
    end
    % normalize
    p_u=p_u/sum_p;  

    n=1;
    % compute the weight value of each pixel in the candidate window
    for i=rmin:rmax
        for j=cmin:cmax
            if (i>=1 & i<=height & j>=1 & j<=width)
                if textureMat(i,j)~=0
                    R=floor(Image(i,j,1)/(256/redBins))+1;
                    G=floor(Image(i,j,2)/(256/greenBins))+1;
                    B=floor(Image(i,j,3)/(256/blueBins))+1;
                    T=textureMat(i,j);
                    u=(R-1)*greenBins*blueBins*5+(G-1)*blueBins*5+(B-1)*5+T;
                    x(1,n)=i;
                    x(2,n)=j;
                    w_i(n)=sqrt(q_u(u)/p_u(u));
                    n=n+1;
                end
            end
        end
    end

    center_r=(x*w_i'/sum(w_i))';               % new center  
    MS=sqrt(sum((center_r-center_old).^2));    % normal of mean shift vector 
    iterations=iterations+1;                   % count iteration number
    
   % Does mean shift algorithm converge?    
    if (MS<minDist | iterations>=maxIterNum)
        break;
    end
    
    center_old=center_r;                        % save tracking result        
    center=floor(center_r);                    

    % reInitialize the candidate winodw
    rmin=center(1)-w_halfsize(1)-incre;        
    rmax=center(1)+w_halfsize(1)+incre;
    cmin=center(2)-w_halfsize(2)-incre;
    cmax=center(2)+w_halfsize(2)+incre;
    
    histo=zeros(redBins,greenBins,blueBins,5);  % reInitialize the target histogram
    sum_p=0;         
end