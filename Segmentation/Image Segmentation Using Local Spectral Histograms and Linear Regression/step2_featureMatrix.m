% Step 2: Compute the feature matrix. The plot shows the singular values. 

ws=9; % integration scale. The actual value is 2*ws+1. 

sh_mx=zeros(bn*bb,N1,N2,'single');

bs=1:bn*bb;
for i=1:N1
    for j=1:N2
        wtl=[max(i-ws,1),max(j-ws,1)];
        wbr=[min(i+ws,N1),min(j+ws,N2)];
        sz=(wbr(1)-wtl(1)+1)*(wbr(2)-wtl(2)+1);
        if wtl(1)==1 && wtl(2)==1
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2));
        elseif wtl(1)==1 && wtl(2)~=1
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2))-HImap2(bs,wbr(1),wtl(2)-1);
        elseif wtl(1)~=1 && wtl(2)==1
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2))-HImap2(bs,wtl(1)-1,wbr(2));
        else
            sh_mx(bs,i,j)=HImap2(bs,wbr(1),wbr(2))+HImap2(bs,wtl(1)-1,wtl(2)-1)...
                -HImap2(bs,wbr(1),wtl(2)-1)-HImap2(bs,wtl(1)-1,wbr(2));
        end
        sh_mx(bs,i,j)=sh_mx(bs,i,j)/sz;
    end
end

Y=reshape(sh_mx,bb*bn,N1*N2);
S=Y*Y';
[v,d]=eig(S);
k=flipud(sqrt(diag(abs(d))));
plot(k,'*')

