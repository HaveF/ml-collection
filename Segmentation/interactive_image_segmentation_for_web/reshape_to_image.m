%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function img = reshape_to_image(X, width, height, flag_assignment)

% each cloumn is a data point

% X:  cooresponding to the image data row by row

% width:            the width
% height:           the height
% flag_assignment:    = 1;                    assign the values to  the 3 * 3 neighbors
%                             = 0 (otherwise);  NOT assign the values to  the  neighbors  


[dim, N] = size(X);

%size_width = floor(sqrt(dim));  % we hope it is a square

size_width = 3;
half_width = 1;

if nargin < 4 
    size_width = 3;
    half_width = floor(size_width / 2);
end

if flag_assignment == 1 
    size_width = 3;                     % assign the regressed 
    half_width = 1;
else
    size_width = 1;
    half_width = 0;
end

img = zeros(height, width);

img_accessed = zeros(height, width);


p = 1;
for i = 1: height
        ystart = i - half_width;
        yend   = i + half_width;
    
        if ystart < 1
            ystart = 1;
            yend = ystart + size_width - 1;
        end
    
        if(yend > height)
            yend   = height;
            ystart = yend - size_width + 1;
        end
     
    for j = 1: width
        
        xstart = j - half_width;
        xend   = j + half_width;
    
    if xstart < 1
        xstart = 1;
        xend = xstart + size_width - 1;
    end
    
    if(xend > width)
        xend   = width;
        xstart = xend - size_width + 1;
    end
        
        T = X(p);                                    % using for scalar mapping   
        T1 = T .* ones(size_width, size_width);      % using for scalar mapping        
        
        img(ystart:yend, xstart:xend) = img(ystart:yend, xstart:xend) + T1;
        
        img_accessed(ystart:yend, xstart:xend) = img_accessed(ystart:yend, xstart:xend) + ones(size_width, size_width);
     
        p = p + 1;
    
    end
    
end

img = img ./ img_accessed;    % calculate the mean by accessment times








