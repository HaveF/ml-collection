% Step 1: Compute the integral histograms

clear
close all

bn=7; % the number of filters

I=imread('6.jpg');
IMG =I;
%I = rgb2gray(I);
%I=imread('fruit.jpg');
LAB = colorspace('LAB<-RGB',double(I));
%I=LAB;
[N1,N2,N3]=size(I);
I=single(I);

Ig=zeros(N1,N2,bn,'single');

% the filterbank (adjust if necessary)
%%  color info 3 %%%%
Ig(:,:,1:3) = LAB;
%% log 2 %%
%Ig(:,:,12:14) = I;
h=fspecial('log',[5,5],.5);
Ig(:,:,end+1) = imfilter(I(:,:,1),h,'symmetric');
h=fspecial('log',[7,7],1);
Ig(:,:,end+1) = imfilter(I(:,:,1),h,'symmetric');
%% position 2 %%
XX = (1:N1)';
XXrep = repmat(XX,1,N2);
Ig(:,:,end+1) = XXrep;
XX = (1:N2);
XXrep = repmat(XX,N1,1);
Ig(:,:,end+1) = XXrep;
% 
% %% gabor 4 %%
% h=gabor_fn(3,pi/2);
% Ig(:,:,end+1) = imfilter(I(:,:,1),h,'symmetric');
% h=gabor_fn(3,0);
% Ig(:,:,end+1) = imfilter(I(:,:,1),h,'symmetric');
% h=gabor_fn(3,pi/4);
% Ig(:,:,end+1) = imfilter(I(:,:,1),h,'symmetric');
% h=gabor_fn(3,-pi/4);
% Ig(:,:,end+1) = imfilter(I(:,:,1),h,'symmetric');
% h=gabor_fn(3,pi/2);
% Ig(:,:,end+1) = imfilter(I(:,:,2),h,'symmetric');
% h=gabor_fn(3,0);
% Ig(:,:,end+1) = imfilter(I(:,:,2),h,'symmetric');
% h=gabor_fn(3,pi/4);
% Ig(:,:,end+1) = imfilter(I(:,:,2),h,'symmetric');
% h=gabor_fn(3,-pi/4);
% Ig(:,:,end+1) = imfilter(I(:,:,2),h,'symmetric');
% h=gabor_fn(3,pi/2);
% Ig(:,:,end+1) = imfilter(I(:,:,3),h,'symmetric');
% h=gabor_fn(3,0);
% Ig(:,:,end+1) = imfilter(I(:,:,3),h,'symmetric');
% h=gabor_fn(3,pi/4);
% Ig(:,:,end+1) = imfilter(I(:,:,3),h,'symmetric');
% h=gabor_fn(3,-pi/4);
% Ig(:,:,end+1) = imfilter(I(:,:,3),h,'symmetric');


bb=16; % bin number

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


