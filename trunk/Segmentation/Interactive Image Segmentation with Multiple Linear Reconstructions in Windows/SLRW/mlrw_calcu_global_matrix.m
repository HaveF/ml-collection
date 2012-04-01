
%    This code is written by Shiming Xiang
% 
%    Address:    National Laboratory of Pattern Recognition (NLPR), 
%                      Institute of Automation, Chinese Academy of Sciences, Beijing, 100190, China
%    Email:        smxiang@gmail.com

%    This software can be used freely for research purposes.
%    Published reports of research  using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%
%    Shiming Xiang, Chunhong Pan, Feiping Nie, Changshui Zhang: 
%    Interactive Image Segmentation With Multiple Linear Reconstructions in Windows. 
%    IEEE Transactions on Multimedia, Vol. 13, No.2, pp. 342-352, 2011.

%    This code has been compiled and tested using matlab 6.5  and matlab     7.0

%    version 1.0, 
%    updated by Shiming Xiang,  Oct. 19, 2011  

%==========================================================================


%==========================================================================
% ABOUT COMPUTATION TIME:
%
% With pure matlab codes, finishing all the computations may cost hundreds  of seconds. 
% For  example,  Finishing all of the pixels related to sentence  "M(ind, ind) = M(ind, ind) + ipl" in this function
% may cost hundreds of seconds on PC. However, this step of adding local matrix
% to global matrix can be easily implemented with C language
% (For example, it  may cost less than 1 sencod  for images with 400 * 400 pixels on most PCs.).   
%==========================================================================



function M = mlrw_calcu_global_matrix(X, neighborhood, lambda)

% X:                         each column is a data point
% neighborhood:       each column records the indices of neighbors, including itself 
% lambda:                the regularization parameter in each neighborhoods


[K, N] = size(neighborhood);   


%caluclate  the reconstruction matrix;
tol = 0.0001;
if nargin >  2
    tol = lambda;
end

M = sparse(1:N, 1:N, zeros(1,N), N, N, 4*K*N); 

for ii=1:N
   
    thisx = X(:, neighborhood(:, ii));     % including itself, each columnn is a data point
        
    ipl = zeros(K, K);
    
    
    for j = 1: K
        index = 1:K;
        index(j) = [];
        z = thisx(:, index) - repmat(thisx(:,j), 1, K-1);   % shift ith pt to origin
    
        %% = =  = = =  = =  = = =  = =  = = =  = =  = = =  = =  = = =  = =
        C = z'*z;                                           
        C = C + eye(K-1,K-1) * lambda;     

        W = C \ ones(K-1,1);                       
                
        W = W / sum(W);                             
        
        %construct the alignment matrix
        TT = zeros(K, K);
        TT(j, j) = 1;
        TT(j, index) = -W';
        TT(index, j) = -W;
        TT(index, index) = W * W';

        ipl = ipl + TT;
    end
    
    ind = neighborhood(:, ii);
    
   M(ind, ind) = M(ind, ind) + ipl;         % Note that this step is very time-consuming. 
                                                           % For example,  with 400 * 400 pixel size image, finishing this step of
                                                           % "M(ind, ind) = M(ind, ind) + ipl"  may cost hundreds of seconds on PC. 
                                                           
                                                           % We can use C  codes to  replace this  step of  assembing; 
                                                           % Note that the computation  time reported  in the paper  is evaluated with C codes  for this step.
                                                           % 
                                                           % Note again that we can also  use C codes to replace all of the above matlab sentences.
                                                           % (For example, it  may cost less than 1 sencod  for images with 400 * 400 pixels on most PCs.)   
   
end;
