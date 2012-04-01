%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function X = cvrt_img_to_matrix(Img)

% Img:    The image to be segmented

% X  :    the pixel information, in decimal form

%convert it to a decimal matrix, each column is a data points

% scan the image row by row, .......

[h, w] = size(Img(:,:,1));
size1 = h * w;

R1 = Img(:,:,1);
G1 = Img(:,:,2);
B1 = Img(:,:,3);

R1 = R1';
G1 = G1';
B1 = B1';

R1 = reshape(R1, size1, 1);
G1 = reshape(G1, size1, 1);
B1 = reshape(B1, size1, 1);

X = [R1,G1,B1]';

X = double(X);

s = 1 / 255.0;

X = s .* double(X);











