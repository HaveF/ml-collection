function [textureMat,textureVarMat]=LBPTexture(X,P,R,T)
% **************************************
% function: calculate lbp texture image of X
% parameters:
%       X: image
%       P,R: pattern of lbp texture 
%       T, threshold
% returns: 
%       textureMat: lbp texture of X
%       VAR: variance of texture    
%
% Note: we refer to  lbp matlab code of University of Oulu about the
% implementation of lbp operator
% **************************************


patternMap = [1 2 2 3 2  10 3 4 2 10  10 10 3 10 4  5 2 10 10 10  10 10 10 10 3  10 10 10 4 10  5 6 ...
	2 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 3 10  10 10 10 10 10  10 4 10 10 10  5 10 6 7 ...
	2  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 ... 
	3 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 4 10 10  10 10 10 10 10  5 10 10 10 6  10 7 8 ...
	2 3  10 4 10 10 10  5 10 10 10 10  10 10 10 6 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 7 ...
	10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 10 10  10 8 ...
	3 4 10  5 10 10 10 6  10 10 10 10 10  10 10 7 10 10  10 10 10 10 10  10 10 10 10 10  10 10 10 8 ...
	4  5 10 6 10 10  10 7 10 10 10  10 10 10 10 8  5 6 10 7 10  10 10 8 6 7  10 8 7 8 8 9];

w1 = (1/sqrt(2))^2;
w2 = (1-1/sqrt(2))*(1/sqrt(2));

D = size(X);
sy = D(1);
sx = D(2);

Xi = zeros(sy+2,sx+2);
Xi(2:sy+1,2:sx+1) = X;

Xi1 = zeros(sy+2,sx+2);Xi2 = zeros(sy+2,sx+2);Xi3 = zeros(sy+2,sx+2);Xi4= zeros(sy+2,sx+2);
Xi5 = zeros(sy+2,sx+2);Xi6 = zeros(sy+2,sx+2);Xi7 = zeros(sy+2,sx+2);Xi8 = zeros(sy+2,sx+2);
p1 = zeros(sy+2,sx+2);p2 = zeros(sy+2,sx+2);p3 = zeros(sy+2,sx+2);
p4 = zeros(sy+2,sx+2);p5 = zeros(sy+2,sx+2);p6 = zeros(sy+2,sx+2);
p7 = zeros(sy+2,sx+2);p8 = zeros(sy+2,sx+2);p9 = zeros(sy+2,sx+2);


p1(3:sy+2,3:sx+2) = X ;
p2(3:sy+2,2:sx+1) = w2*double(X) ;
p3(3:sy+2,1:sx) = X ;

p4(2:sy+1,3:sx+2) = w2*double(X) ;
p5(2:sy+1,2:sx+1) = (1-1/sqrt(2))^2*double(X) ;
p6(2:sy+1,1:sx) = w2*double(X) ;

p7(1:sy,3:sx+2) = X ;
p8(1:sy,2:sx+1) = w2*double(X) ;
p9(1:sy,1:sx) = X ;

Xi1 = w1*p1+ p2+p4 + p5 + 0.000001; %Xi1 to the right and down from X
Xi2(3:sy+2,2:sx+1) = X;
Xi3 = w1*p3 + p2 + p6 + p5 + 0.000001;
Xi4(2:sy+1,1:sx) = X;
Xi5 = w1*p9 + p8 + p6 + p5 + 0.000001;
Xi6(1:sy,2:sx+1) = X;
Xi7 = w1*p7 + p8 + p4 + p5 + 0.000001;
Xi8(2:sy+1,3:sx+2) = X;

% calculate mean
meanXi=(Xi1+Xi2+Xi3+Xi4+Xi5+Xi6+Xi7+Xi8)/8;
% calculate variance
VAR=sqrt(((Xi1-meanXi).^2+(Xi2-meanXi).^2+(Xi3-meanXi).^2+(Xi4-meanXi).^2+(Xi5-meanXi).^2+(Xi6-meanXi).^2+(Xi7-meanXi).^2+(Xi8-meanXi).^2))/8;

% calculate lbp texture
Xi= (Xi4>=Xi+T)+2*(Xi5>=Xi+T)+4*(Xi6>=Xi+T)+8*(Xi7>=Xi+T)+16*(Xi8>=Xi+T)+32*(Xi1>=Xi+T)+64*(Xi2>=Xi+T)+128*(Xi3>=Xi+T);
X=Xi(3:sy,3:sx);
M=size(X,1);
N=size(X,2);

for i=1:M
    for j=1:N
        textureMat(i,j)=patternMap(X(i,j)+1);        
    end
end
textureVarMat=VAR(3:sy,3:sx);  
avg=sum(textureVarMat(:))/(N*M);
textureVarMat=round(255*textureVarMat/max(textureVarMat(:)));

return;