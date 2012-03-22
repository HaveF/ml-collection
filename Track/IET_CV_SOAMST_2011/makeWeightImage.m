function makeWeightImage(image,x0,y0,Cov,p_u,q_u,frameNo);

k=1.2;
[U S V]=svd(Cov);
S(1,1)=(sqrt(S(1,1)*k))^2;
S(2,2)=(sqrt(S(2,2)*k))^2;
Cov=U*S*V;

Ind=moment2Ellipse(Cov);
Ind(:,1)=Ind(:,1)+x0;
Ind(:,2)=Ind(:,2)+y0;

p_u=rgbPDF_SOAMST(image,x0,y0,Ind(:,1),Ind(:,2),16,16,16);

N=size(Ind,1);

w_u=sqrt(q_u./(p_u+eps));

for i=1:N
    row=Ind(i,2);
    col=Ind(i,1);

    R=floor(image(row,col,1)/(256/16))+1;
    G=floor(image(row,col,2)/(256/16))+1;
    B=floor(image(row,col,3)/(256/16))+1;
    u=(R-1)*16*16+(G-1)*16+B;
    
    weiImg(row,col)=w_u(u);
end