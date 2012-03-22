function image=makeTrackingImage(image,centerX,centerY,Ind,colorType)

H=size(image,1);
W=size(image,2);

N=size(Ind,1);

for i=1:N
    col=Ind(i,1)+centerX;
    row=Ind(i,2)+centerY;
    if (row>=1 & row<=H & col>=1 & col<=W)
        switch colorType
            case 1
                image(row,col,:)=255;        
            case 2
                image(row,col,:)=0;
            case 3
                image(row,col,1)=255;
                image(row,col,2)=0;
                image(row,col,3)=0;
            case 4
                image(row,col,1)=0;
                image(row,col,2)=255;
                image(row,col,3)=0;
            case 5
                image(row,col,1)=0;
                image(row,col,2)=0;
                image(row,col,3)=255;
            case 6
                image(row,col,1)=255;
                image(row,col,2)=255;
                image(row,col,3)=0;
            case 7
                image(row,col,1)=0;
                image(row,col,2)=255;
                image(row,col,3)=255;
            case 8
                image(row,col,1)=255;
                image(row,col,2)=0;
                image(row,col,3)=255;
        end      
    end
end