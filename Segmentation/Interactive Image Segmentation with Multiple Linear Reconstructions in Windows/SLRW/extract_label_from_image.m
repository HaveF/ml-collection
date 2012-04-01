%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 

function  Inds = extract_label_from_image(Img, mask_color)
% Img:          [H,W,3]
% mask_color:   the mask colors,each row is a mask color with RGB order; size = number_of_vision_objects * 3; 

% return
% Inds:          the extracted indeces of masked colors, one column corresponds to one mask 

%here, we assume that the label information is drawn in the image
%      We also assume that the pixels are indexed  in a row-scanning way.

RR = Img(:,:,1);
GG = Img(:,:,2);
BB = Img(:,:,3);

RR = RR';       % Note that Matlab indentify the matrix elements column  by column
GG = GG';
BB = BB';

number_mask = size(mask_color, 1);
Inds = cell(number_mask, 1);

for i= 1: number_mask
    Ind_temp = find(RR == mask_color(i, 1) & GG == mask_color(i, 2) & BB == mask_color(i, 3) );
    Inds{i, 1}  = Ind_temp;    
end


return;
