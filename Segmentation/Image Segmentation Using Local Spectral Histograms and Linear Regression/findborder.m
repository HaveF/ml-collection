function [result connect_matrix]= findborder(A)
    %A= rand(50,50);
    [N1,N2] = size(A);
    maxregionno = max(max(A));
    connect_matrix = zeros(maxregionno,maxregionno);
    result = zeros(N1,N2);
    result = uint8(result);
    result(1:2,:) = 255;
    result(N1-1:N1,:) = 255;
    result(:,1:2) = 255;
    result(:,N2-1:N2) = 255;
    for r=2:(N1-1)
        for c = 2:(N2-1)
%             if A(r,c)~=A(r-1,c-1)|A(r,c)~=A(r-1,c)|A(r,c)~=A(r-1,c+1)...
%                     |A(r,c)~=A(r+1,c)|A(r,c)~=A(r+1,c+1)|A(r,c)~=A(r,c+1)...
%                     |A(r,c)~=A(r+1,c-1)|A(r,c)~=A(r-1,c+1)
%                 result(r,c) = 255;
%             end

            for step1=-1:1
                for step2 =-1:1
                    if A(r,c)~=A(r+step1,c+step2)
                        result(r,c) = 255;
                        connect_matrix(A(r,c),A(r+step1,c+step2))=1;
                    end
                end
            end
        end
    end
end