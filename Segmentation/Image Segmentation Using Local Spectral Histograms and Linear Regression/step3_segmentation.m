% Step 3: produce segmentation

segn=5; % the number of segments 
edgeThr=.5; % the threshold for removing edge features 
            % may need to adjust if the filterbank is changed

dimn=segn;
v1=fliplr(v);
U1=v1(:,1:dimn);

Y2=Y;
Y1=(Y2'*U1)';

intreg=zeros(N1,N2);
for i=1+ws:N1-ws
    for j=1+ws:N2-ws
        up=U1'*sh_mx(:,i-ws,j);
        bt=U1'*sh_mx(:,i+ws,j);
        lf=U1'*sh_mx(:,i,j-ws);
        rt=U1'*sh_mx(:,i,j+ws);
%         tmp=sum((up-bt).^2./(up+bt+eps))+sum((lf-rt).^2./(lf+rt+eps));
        tmp=sqrt(sum((up-bt).^2)+sum((lf-rt).^2));
        if tmp<edgeThr
            intreg(i,j)=1;
        end
    end
end 

imshow(intreg,[])
title('removing features on edges')
%pause
idx=find(intreg==1);
len=length(idx);
Mx=Y1(:,idx);

tmplt=zeros(dimn,segn,'single');
L=sum(Mx.^2);
[y rn]=max(L);
tmplt(:,1)=Mx(:,rn); 
n=1; 

seedmap=zeros(N1,N2);
seedmap(idx(rn))=1;

tn=n+1;
CY=repmat(tmplt(:,n),1,len);
ccos=sqrt(sum((Mx-CY).^2));
[m id]=max(ccos);
tmplt(:,tn)=Mx(:,id);
seedmap(idx(id))=1;
    
while tn<segn
    tmp=zeros(tn,len);
    for i=1:tn
        CY=repmat(tmplt(:,i),1,len);
        tmp(i,:)=sqrt(sum((Mx-CY).^2));
    end
    tn=tn+1;
    ccos=min(tmp); 
    [m id]=max(ccos);
    tmplt(:,tn)=Mx(:,id);
    seedmap(idx(id))=1;
end
  
% seedmap=imdilate(seedmap,ones(5));
% imshow(seedmap,[])
% pause

% Kmeans clustering

cenInt=tmplt;
ccos=zeros(segn,len,'single');
flag=1;

while flag==1    
for i=1:segn
    CY=repmat(cenInt(:,i),1,len);
    ccos(i,:)=sqrt(sum((Mx-CY).^2));
end

[M clab]=min(ccos);
NcenInt=zeros(dimn,segn,'single');

for i=1:segn
    tmind=find(clab==i);
    tmp=Mx(:,tmind);
    NcenInt(:,i)=sum(tmp,2)./length(tmind);
end

if NcenInt==cenInt
    flag=0;
else
    cenInt=NcenInt;
end
end

[M clab]=min(ccos);
kmres=zeros(N1,N2);
kmres(idx)=clab+1;
% figure,imshow(kmres,[])
% pause

B=(NcenInt'*NcenInt)^(-1)*NcenInt'*Y1;

[m,slab]=max(B);

finres=reshape(slab,N1,N2);

figure,imshow(finres,[],'border','tight')


