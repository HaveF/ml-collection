%    version 1.0, 
%    rewritten and updated by Shiming Xiang,  Oct. 6, 2010  
function nb = calcu_only_spatial_neighorhood(width, height, win_width, win_height)

% width:      image with;
% height:     image height
% win_width:  the width  of the window for calculating the neighbors
% win_height: the height of the window for calculating the neighbors

% nb:         the neighbor indces, recored in each column. In each column, the first is itself, and then its their neighbors 

N  = width * height;
K = win_width * win_height;  % the number of the neighbors 

nb = zeros(K, N);

index = ones(K, 1);

w_half = floor(win_width / 2);
h_half = floor(win_height / 2);

upper_height = height - win_height + 1;
upper_width  = width  - win_width + 1;

%Find the neighborhood

pp= 1;
for i = 1: height
    start_h = i - h_half;
    
    if start_h < 1
        start_h = 1;
    end
    if start_h > upper_height
        start_h = upper_height;
    end
            
    for j = 1: width
        
        start_w = j - w_half;
    
        if start_w < 1
            start_w = 1;
        end
    
        if start_w > upper_width
            start_w = upper_width;
        end 
        
        index(1) = (i - 1) * width + j;   % the center point, thus including itself
        p = 2;
        
        
        for ii = start_h : start_h + win_height - 1
            for jj = start_w : start_w + win_width - 1
                if ii ~= i | jj ~= j
                    index(p) = (ii - 1) * width + jj;
                     p = p+1;
                end
            end
        end
    
        nb(:, pp) = index;
        pp = pp + 1;
        
    end 
  
end

