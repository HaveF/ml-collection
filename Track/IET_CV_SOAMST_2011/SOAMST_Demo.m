% =========================================================================
% Filename: SOAMST_Demo.m                                                 *
% Function: Implement the algorithm described in "Scale and Orientation   *
%           Adaptive Mean Shift Tracking.",                               *
%           Accepted by IET Computer Vision, 2011.                        *
% Author:   Jifeng Ning, Lei Zhang, David Zhang, Chengke Wu              *
% email:    jf_ning@sina.com,cslzhang@comp.polyu.edu.hk                   *        
% Date  :   March,8,2011                                                  *
% All Copyright Reserved.                                                 *
%                                                                         *
% Address   State Lab on Intergrated Server Network, Xidian Univeristy     *
%           Department of Computing, The Hong Kong Polytechnic University  *
%           College of Information Engineering, Northwest A&F University  *                       
%                                                                         *
% =========================================================================

clc;
clear all;

aviName='ellipse_slow';
avi = aviread(aviName);                       % load avi file
[dontneed numberofframes] = size(avi);        % numberofframes, number of frame
Frames={avi.cdata};                           % frames, video data    
frame0=Frames{1};                             % first frame    

redBins=16;greenBins=16;blueBins=16;          % default color quant sheme (fixed in this demo)

startFrame=1;                                 %  
endFrame=numberofframes;                      % 
incre=input('Increased search size (5-15 is good): ');
omiga=input('Corrected area pararmeter (1-2 is good): ');

figure(1),imMask=roipoly(frame0);             %
[ry,rx]=find(imMask);                         % 
m0=floor(mean([rx,ry]))';                     % 
V0=cov([rx,ry]);                              % 
V0=correctCov(V0,size(rx,1));                 % 
q_u=rgbPDF_SOAMST(double(frame0),m0(1),m0(2),rx,ry,redBins,greenBins,blueBins);  % calculate the target model

Cov=V0;                                       % 
V0_searchRegion=enlargeCov(Cov,incre);        % 
Cov_searchRegion=V0_searchRegion;             % 

imgHeight=size(frame0,1);                     % 
imgWidth=size(frame0,2);                      % 

x_0=m0(1);y_0=m0(2);                          % 

saveResults=1;                                % save results
if saveResults==1
    avi=avifile('E:\SOAMST7.avi');
end

for i=startFrame:endFrame   
    [Ind_searchRegion,Ind_searchRegion_EpanchnikovKernel]=moment2Ellipse(Cov_searchRegion);    
    [x_0,y_0,CovNew,BhattCoff,iterations,p_u]=rgbTracking_SOAMST(double(Frames{i}),imgWidth,imgHeight,...
    q_u,x_0,y_0,...
    Ind_searchRegion,Ind_searchRegion_EpanchnikovKernel,...
    omiga,redBins,greenBins,blueBins,i);

    iterNum(i)=iterations; 
    
    figure(1);

    Ind_targetRegion_AVI=moment2Ellipse(CovNew);                                                
    Ind_searchRegion_AVI=moment2Ellipse(Cov_searchRegion);                                      
    
    Ind_targetRegionBoundary_AVI=findRegionBoundary_AVI(Ind_targetRegion_AVI);                  
    Ind_searchRegionBoundary_AVI=findRegionBoundary_AVI(Ind_searchRegion_AVI);                   
    
    frame=uint8(makeTrackingImage(double(Frames{i}),x_0,y_0,Ind_targetRegionBoundary_AVI,3));    
    frame=uint8(makeTrackingImage(frame,x_0,y_0,Ind_searchRegionBoundary_AVI,4));               
    
    imshow(frame);title(['SOAMST ',num2str(i),'/',num2str(numberofframes)]),drawnow; % display tracking results
    
    Cov_searchRegion=enlargeCov(CovNew,incre);     %
    if saveResults==1
        avi=addframe(avi,frame);
    end
end

if saveResults
    avi=close(avi);
end