% Step 1: Compute the integral histograms

clear
close all

bn=3; % the number of filters

I=imread('texturetest_resized.jpg');
%I=imread('fruit.jpg');

[N1,N2]=size(I);
I=single(I);

Ig=zeros(N1,N2,bn,'single');

% the filterbank (adjust if necessary)
Ig(:,:,1) = I;
h=fspecial('log',[5,5],.5);
Ig(:,:,2) = imfilter(I,h,'symmetric');
h=fspecial('log',[7,7],1);
Ig(:,:,3) = imfilter(I,h,'symmetric');

% h=gabor_fn(3,pi/2);
% Ig(:,:,4) = imfilter(I,h,'symmetric');
% h=gabor_fn(3,0);
% Ig(:,:,5) = imfilter(I,h,'symmetric');
% h=gabor_fn(3,pi/4);
% Ig(:,:,6) = imfilter(I,h,'symmetric');
% h=gabor_fn(3,-pi/4);
% Ig(:,:,7) = imfilter(I,h,'symmetric');

bb=11; % bin number

HImap2=zeros(bn*bb,N1,N2,'single');

for b=1:bn
    Ic=single(Ig(:,:,b));
    mx=max(Ic(:));
    mn=min(Ic(:));
    U=(mx-mn)/bb;
    binc=mn+(1:bb)*U-U/2;


    difc=abs(Ic(1,1)-binc);
    [d,id]=min(difc);
    tmp=zeros(1,bb,'single');
    tmp(id)=1;
    HImap2((bb*b-bb+1):bb*b,1,1)=tmp';

    for j=2:N2
        difc=abs(Ic(1,j)-binc);
        [d,id]=min(difc);
        tmp(id)=tmp(id)+1;
        HImap2((bb*b-bb+1):bb*b,1,j)=tmp';
    end

    for i=2:N1
        difc=abs(Ic(i,1)-binc);
        [d,id]=min(difc);
        HImap2(bb*b-bb+id,i,1)=HImap2(bb*b-bb+id,i-1,1)+1;
        tmp=zeros(1,bb,'single');
        tmp(id)=1;
        for j=2:N2
            difc=abs(Ic(i,j)-binc);
            [d,id]=min(difc);
            tmp(id)=tmp(id)+1;
            HImap2((bb*b-bb+1):bb*b,i,j)=HImap2((bb*b-bb+1):bb*b,i-1,j)+tmp';
        end
    end
end



