%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 



function MM = post_process_image_for_two_class_regression(X_test_new, width, height, flag_assignment)

%  X_test_new:          the regressed results, a row vector
%  [width, height]:       the width and height of the image
% flag_assignment:       = 1;                    assign the values to  the 3 * 3 neighbors
%                                = 0 (otherwise);  NOT assign the values to  the  neighbors  


%  return: 
%  the segmented results

if nargin < 4 
    flag_assignment = 1;      % do the neighbor assignment
end

% now we need to add them together
N = width * height;
MM = reshape_to_image(X_test_new, width, height, flag_assignment);
Y_label = reshape(MM, N, 1);


%=========================================================================
% construct the threshold
threshold = 0.0;

%=========================================================================
% using different graysacles to show the background and foreground objects

YY = zeros(N, 1);
index = find(Y_label > 0.0 );
YY(index) = 255;

MM = reshape( YY, height, width);

MM = uint8(MM);

% Here I suggest to add the following sentence, Apr.21, 2010, by Shiming Xiang
MM = medfilt2(MM, [3 3]); 
