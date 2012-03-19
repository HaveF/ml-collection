function N1=SmoothContourNormals2D(N1,Lines1,n)
if(nargin<3), n=10; end
for k=1:n
    Nupdate=zeros(size(N1));
    Nupdate(Lines1(:,1),:)=N1(Lines1(:,2),:);
    Nupdate(Lines1(:,2),:)=N1(Lines1(:,1),:);
    N1=N1+Nupdate/4;
    L=sqrt(N1(:,1).^2+N1(:,2).^2);
    N1(:,1)=N1(:,1)./L;
    N1(:,2)=N1(:,2)./L;
end