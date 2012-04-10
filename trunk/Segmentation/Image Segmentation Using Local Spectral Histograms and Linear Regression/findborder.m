function result = findborder(A)
    %A= rand(50,50);
    [N1,N2] = size(A);
    result = zeros(N1,N2);
    result = uint8(result);
    result(1,:) = 255;
    result(N1,:) = 255;
    result(:,1) = 255;
    result(:,N2) = 255;
    for r=2:(N1-1)
        for c = 2:(N2-1)
            if A(r,c)~=A(r-1,c-1)|A(r,c)~=A(r-1,c)|A(r,c)~=A(r-1,c+1)...
                    |A(r,c)~=A(r+1,c)|A(r,c)~=A(r+1,c+1)|A(r,c)~=A(r,c+1)...
                    |A(r,c)~=A(r+1,c-1)|A(r,c)~=A(r-1,c+1)
                result(r,c) = 255;
            end
        end
    end
end