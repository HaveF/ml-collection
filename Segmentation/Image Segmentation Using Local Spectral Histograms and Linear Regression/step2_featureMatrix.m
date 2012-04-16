% Step 2: Compute the feature matrix. The plot shows the singular values. 

ws=9; % integration scale. The actual value is 2*ws+1. 

sh_mx=zeros(bn*bb,N1,N2,'single');
bs=1:bn*bb;
for i=1:N1
    for j=1:N2
        wtl=[max(i-ws,1),max(j-ws,1)];   %top left point
        wbr=[min(i+ws,N1),min(j+ws,N2)]; %bottom right point
        sz=(wbr(1)-wtl(1)+1)*(wbr(2)-wtl(2)+1); %area,windows size
        if wtl(1)==1 && wtl(2)==1
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2));      %margin point
        elseif wtl(1)==1 && wtl(2)~=1
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2))-HImap2(bs,wbr(1),wtl(2)-1); %margin point
        elseif wtl(1)~=1 && wtl(2)==1
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2))-HImap2(bs,wtl(1)-1,wbr(2));%margin point
        else
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2))+HImap2(bs,wtl(1)-1,wtl(2)-1)...
                -HImap2(bs,wbr(1),wtl(2)-1)-HImap2(bs,wtl(1)-1,wbr(2));   %fig. 1
        end
        sh_mx(bs,i,j)=sh_mx(bs,i,j)/sz;
    end
end

Y=reshape(sh_mx,bb*bn,N1*N2);
% Y1 = f_cal_standardized_feature(Y')';
% Y=Y1;
S=Y*Y';
[v,d]=eig(S);
k=flipud(sqrt(diag(abs(d))));
plot(k,'*')

