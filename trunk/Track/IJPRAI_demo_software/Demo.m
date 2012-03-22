% *************************************************************************
% Filename: Demo.m                                    *
% Function: Implement the algorithm described in "Robust Object Tracking
% using Joint Color-Texture Histogram"
% International Journal of Pattern Recognition and Artifical Intelligence,vol.23,no.7,pp.1245-1263, 2009.
% Author    Jifeng Ning, Lei Zhang, David Zhang, Chengke Wu     
% Email;    cslzhang@comp.polyu.edu.hk
% Date    : Nov. 25, 2009                                                
% Address   College of Information Engineering, Northest A&F University
%           State Key Laboratory of Integrated Service Networks, Xidian University 
%           Department of Computing, The Hongkong Polytechnic University                                        
%                                                                         
%  Copyright (c) 2009 by Jifeng Ning, Lei Zhang, David Zhang, Chengke Wu  
% *************************************************************************

clear all;
clc;

M = aviread('tabletennis.avi');          % Read video file
[dontneed numberofframes] = size(M);     % the number of frames of the video
Frames={M.cdata};              

% select tracking algorithm
% algorithm: 1: mean shift tracking with rgb histogram  
%            2: mean shift tracing with proposed color-texute joint histogram
algorithm=2;        % default select

minDist=0.1;        % the minimum distance between two mean shift iteration
maxIterNum=15;      % maximal iteration number
incre=5;            % enlarge the search region

% color quantification
redBins=8; greenBins=8; blueBins=8;

% lbp threshold
lbpThreshold=8;

% set double buffer
set(gcf,'DoubleBuffer','on');

% set the domain to tracking
startFrm=1;                                  % 
endFrm=numberofframes;                       % 

frame00=Frames{startFrm};                    % get start-frame to select to tracking object
height=size(frame00,1);
width=size(frame00,2);

% select tracking window manually
%[ cmin, cmax, rmin, rmax ] = select( frame00);

%cmin=168; cmax=182; rmin=78; rmax=97;       % for football sequence
cmin=136; cmax=164; rmin=27; rmax=58;        % for table tennis sequence

center(1,1)=floor((rmin+rmax+1)/2);          % the center of window
center(1,2)=floor((cmin+cmax+1)/2);          %                   
            
w_halfsize(1) = round(abs(rmax - rmin)/2);   % half height of window  
w_halfsize(2) = round(abs(cmax - cmin)/2);   % half width of window

% calculate the target model
if algorithm==1      % target model with rgb histogram
    q_u=rgbPDF(double(frame00),center,w_halfsize,redBins,greenBins,blueBins);
elseif algorithm==2  % target model with rgb+lbp histogram
    q_u=flbp81_rgb_PDF(double(frame00),center,w_halfsize,lbpThreshold,redBins,greenBins,1);
end

for i = startFrm:1:endFrm                    % numberofframes        
    if i==startFrm                           % 
        framei=frame00;
    elseif i>1  % 
        % get the current frame
        framei=Frames{i};
        if algorithm==1                      % stand mean shift tracking algorithm
            center=rgbTracking(double(framei),center,w_halfsize,q_u,redBins,greenBins,blueBins,minDist,maxIterNum,incre);        % means shiftËã·¨½øÐÐ¸ú×Ù
        elseif algorithm==2                  % mean shift tracking with the proposed target representation method
            [center,iterations]=flbp81_rgb_Tracking(double(framei),center,w_halfsize,q_u,redBins,greenBins,1,minDist,maxIterNum,incre,lbpThreshold); 
        end
    end
    
    % window corresponding to tracking result
    rmin=center(1)-w_halfsize(1);             
    rmax=center(1)+w_halfsize(1);            
    cmin=center(2)-w_halfsize(2);           
    cmax=center(2)+w_halfsize(2);           
    % 
    [rmin,rmax,cmin,cmax]=normWindow(rmin,rmax,cmin,cmax,height,width);
    
    % making tracking result
    trackim=framei;                            
    for r= rmin:rmax
        trackim(r, cmin-1:cmin,:) = 255;       
        trackim(r, cmax:cmax+1,:) = 255;        
    end
    for c= cmin:cmax
        trackim(rmin-1:rmin, c,:) = 255;
        trackim(rmax:rmax+1, c,:) = 255;        
    end   
  
    %¡¡display tracking result   
    image(trackim);title([num2str(i),'/',num2str(numberofframes)]);drawnow;
end