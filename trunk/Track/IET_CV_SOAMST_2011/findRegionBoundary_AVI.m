function B=findRegionBoundary_AVI(Ind)

max_Ind=max(Ind(:));

delta=max_Ind+5;
Ind=Ind+delta;

tempMat=zeros(2*delta-1,2*delta-1);
N=size(Ind,1);

for i=1:N
    row=Ind(i,2);
    col=Ind(i,1);
    tempMat(row,col)=1;
end

SE=strel('disk',2);
tempMat=tempMat-imerode(tempMat,SE);
[Indy,Indx]=find(tempMat);
Indy=Indy-delta;
Indx=Indx-delta;
B=[Indx,Indy];