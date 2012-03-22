function [Ind,Dist]=moment2Ellipse(Cov)

InvCov=inv(Cov);            % 

[U S V]=svd(Cov);           % 
Len=round(sqrt(S(1,1)))+2;  % 

k=0;
for i=-Len:Len              % 
    for j=-Len:Len
        % 
        t=[i,j]*InvCov*[i,j]';
        if t<=1
            k=k+1;
            Ind(k,1)=i;            
            Ind(k,2)=j;            
            Dist(k)=t;             
        end
    end
end
maxDist=max(Dist(:));
Dist=(maxDist-Dist)/maxDist;       

minDist=min(Dist(find(Dist~=0)));  

ind=find(Dist==0);                  
Dist(ind)=minDist;         