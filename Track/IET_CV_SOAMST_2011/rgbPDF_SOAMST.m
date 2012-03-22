function q_u=rgbPDF_SOAMST(image,cenX,cenY,indX,indY,redBins,greenBins,blueBins)
% function: calculate target model q_u
% parameters:
%       image                        current frame
%       center                       initial position of the target
%       Indx,Indy                    coordinates of the target region         
%       redBins,greenBins,blueBins   quantified parameters
% return:
%       q_u          ¡¡              target model%
% Date: March,7,2011
% ****************************************************************


sum_q=0;
histo=zeros(redBins,greenBins,blueBins);

rb=256/redBins;
gb=256/greenBins;
bb=256/blueBins;

fDist=(indX-cenX).^2+(indY-cenY).^2;
maxD=max(fDist)+1;
fDist=maxD-fDist;          

N=size(indX,1);         

for i=1:N
    R=floor(image(indY(i),indX(i),1)/rb)+1;
    G=floor(image(indY(i),indX(i),2)/gb)+1;
    B=floor(image(indY(i),indX(i),3)/bb)+1;
    histo(R,G,B)=histo(R,G,B)+fDist(i);
end
    
for i=1:redBins
    for j=1:greenBins
        for k=1:blueBins
            index=(i-1)*greenBins*blueBins+(j-1)*blueBins+k;            
            q_u(index)=histo(i,j,k);
            sum_q=sum_q+q_u(index);
        end
    end
end

q_u=q_u/sum_q;   % normalize
return;