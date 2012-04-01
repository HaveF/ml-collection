%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 



function CM = solve_linear_system(X, Y)

% X:     the input matrix, each column is a data point;
% Y:     the destination matrix: each row is a label vector


% return
% CM:    coefficient matrix, each row is a coefficient vector, corresponding to spline function 


LL = construct_coefficient_matrix(X);       % calculate the coefficient matrix of the splines
                                                                % Note that,  the first is   the kernel term,  the second is the linear term

[dim, N] = size(X);
K = size(Y, 2);

YY = [Y; zeros(dim + 1, K)];
CM = LL \ YY;

CM = CM';


return;



