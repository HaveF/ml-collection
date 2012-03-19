function [Err,Udet]=DiscreteJacobian(B)
if(ndims(B)==3)
    [Uxy,Uxx] = gradient(B(:,:,1));
    [Uyy,Uyx] = gradient(B(:,:,2));

    % Loop through all pixel locations
    Udet=(Uxx+1).*(Uyy+1)-Uyx.*Uxy;
else
    [Uxy,Uxx,Uxz] = gradient(B(:,:,:,1));
    [Uyy,Uyx,Uyz] = gradient(B(:,:,:,2));
    [Uzy,Uzx,Uzz] = gradient(B(:,:,:,3));

    % Loop through all pixel locations
    Udet=(Uxx+1).*(Uyy+1).*(Uzz+1) + Uxy.*Uyz.*Uzx + Uxz.*Uyx.*Uzy - ...
         Uxz.*(Uyy+1).*Uzx - Uxy.*Uyx.*(Uzz+1) - (Uxx+1).*Uyz.*Uzy;
end
Err=any(Udet(:)<eps);






