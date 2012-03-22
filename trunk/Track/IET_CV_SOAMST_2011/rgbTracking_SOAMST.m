function [x_0,y_0,Cov,BhattCoff,iterations,p_u]=rgbTracking_SOAMST(image,imgWidth,imgHeight,...
    q_u,centerX,centerY,...
    Ind_candidate,Ind_kernelCoff,...
    omiga,redBins,greenBins,blueBins,frameNo)
% ********************************************************************
% Function:Implement SOAMST tracking method
% Parameters:image: current frame
%      imgWidth,imgHeight:
%      q_u: Target model
%      centerX, centerY: intial center
%      Ind_candidate:    position coordinate
%      Ind_kernelCoff:   coefficient corresponding to each position coordinate
%      omiga:            area parameter
%      redBins,greenBins,blueBins:
%      frameNo:          number of current frame
% return:
%      x_0,y_0: the center of tracking result region
%      Cov: covariance, representing the tracking region
%      BhattCoff: similarity between target candidate and target model
%      iterations: iteration number
%      p_u: target candidate
%
% Author:Jifeng Ning
% e-mail: jf_ning@sina.com
% All Copyright Reserved.
%
% *************************** end *************************************

sum_p=0; %
histo=zeros(redBins,greenBins,blueBins); % feature spaces,quantifized bins are 16-16-16
iterations=0;

x_0=centerX;
y_0=centerY;

x_0r=centerX;
y_0r=centerY;

delta_x=0;
delta_y=0;

N=size(Ind_candidate,1);
while 1
    for n=1:N
        %
        col=Ind_candidate(n,1)+x_0;
        row=Ind_candidate(n,2)+y_0;
        %
        if (row>=1 & row<=imgHeight & col>=1 & col<=imgWidth)
            R=floor(image(row,col,1)/(256/redBins))+1;
            G=floor(image(row,col,2)/(256/greenBins))+1;
            B=floor(image(row,col,3)/(256/blueBins))+1;
            histo(R,G,B)=histo(R,G,B)+Ind_kernelCoff(n);
        end
    end

    for i=1:redBins
        for j=1:greenBins
            for k=1:blueBins
                index=(i-1)*greenBins*blueBins+(j-1)*blueBins+k;
                p_u(index)=histo(i,j,k);
                sum_p=sum_p+p_u(index);
            end
        end
    end

    p_u=p_u/sum_p;

    w_i=zeros(1,N);
    x=zeros(2,N);
    for n=1:N
        col=Ind_candidate(n,1)+x_0;
        row=Ind_candidate(n,2)+y_0;

        if (row>=1 & row<=imgHeight & col>=1 & col<=imgWidth)
            R=floor(image(row,col,1)/(256/redBins))+1;
            G=floor(image(row,col,2)/(256/greenBins))+1;
            B=floor(image(row,col,3)/(256/blueBins))+1;
            u=(R-1)*greenBins*blueBins+(G-1)*blueBins+B;
            x(1,n)=col;
            x(2,n)=row;
            a=sqrt(q_u(u));
            b=sqrt(p_u(u));
            w_i(n)=a/b;
        end
    end

    y_1=x*w_i'/sum(w_i);
    delta_x=y_1(1)-x_0r;  delta_y=y_1(2)-y_0r;
    iterations=iterations+1;

    if (sqrt(delta_x^2+delta_y^2)<0.05 | iterations>15)
        break;
    end

    x_0r=y_1(1);
    y_0r=y_1(2);


    x_0=round(x_0r);
    y_0=round(y_0r);

    histo=zeros(redBins,greenBins,blueBins);
    sum_p=0;
end

BhattCoff=sqrt(q_u)*sqrt(p_u)';

sumW_i=sum(w_i(:));

Dxx=w_i*((x(1,:)-x_0)'.^2)/sumW_i;
Dyy=w_i*((x(2,:)-y_0)'.^2)/sumW_i;
Dxy=w_i*((x(1,:)-x_0)'.*(x(2,:)-y_0)')/sumW_i;

k=(BhattCoff-1)/omiga;
c=exp(k);

area=sumW_i*c;

C=[Dxx,Dxy;Dxy,Dyy];
[U,S,V]=svd(C);
l1=sqrt(S(1,1));
l2=sqrt(S(2,2));

s2=area/(pi*l1*l2);

S=s2*S;
Cov=U*S*V;
return;