% This is an Example of using Shape-Context to create a corresponding
% point model of a set of 10 Mandible 3D-volumes.
%
% It loads 10 datasets and convert them to a point
% description of the outer surface. ShapeContext is used to 
% iteratively register one of the datasets to all other datasets.
% The warping between datasets is kept diffeomorphic, by constraining
% the jacobian of the b-spline transformation grid.
% 
% This example takes about 7 hours with a  Intel(R) 
% Core(TM) i7 CPU 870 @ 2.93GHz with 16,0 GB RAM, Windows 7 and 
% Matlab R2010b 64bit 
% 
% To use the constructed point model to create an ASM/AAM
% run the examples in ActiveModels_version7
%
% Example is written by D.Kroon University of Twente
% March 2011 updated January 2012

% Add all the needed folders and functions to Matlab search path
functionname='Example3D.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'functions'])
addpath([functiondir 'ShapeContext_version2'])
addpath([functiondir 'nonrigid_version25'])
addpath([functiondir 'polygon2voxel_version1j'])
addpath([functiondir 'DistanceFieldCP_version1'])

% This will compile the c-code to fast multi-threaded mex-files. The 
% example also works without the compile c-files but will be slower.
cd([functiondir 'ShapeContext_version2']);
try 
   compile_c_files, 
catch ME
   disp(ME); 
end
cd([functiondir 'polygon2voxel_version1j'])
try 
    mex('polygon2voxel_double.c');
catch ME
   disp(ME); 
end
cd([functiondir 'nonrigid_version25']);
try 
   compile_c_files, 
catch ME
   disp(ME); 
end
cd(functiondir);

% Load one dataset as Master
filename='images3D\segm004.mat';
% Convert the segmented image to contour data
[I1,Faces1,Vertices1,Faces1l,Vertices1l]=GetImageSurface(filename);

% Show the dataset
figure, hold on;
FV.faces=Faces1l; FV.vertices=Vertices1l;
patch(FV,'facecolor',[0 0 1],'edgecolor','none'); camlight; view(3)
plot3([Vertices1(:,1) Vertices1(:,1)]',[Vertices1(:,2) Vertices1(:,2)]',[Vertices1(:,3) Vertices1(:,3)]','r.','markersize',2);
drawnow

for j=1:10,
    close all
    % Load a dataset 
    strj=['000' num2str(j)]; strj=strj(end-2:end);
    filename=['images3D\segm' strj '.mat'];
    disp(['Load Dataset : ' filename]); drawnow('expose');
    % Convert the segmented image to contour data
    [In,Facesn,Verticesn,Facesnl,Verticesnl]=GetImageSurface(filename);

    % Itterative register the contours using shape context
    % regulize the matching with a b-spline grid which is refined from
    % coarse to fine during this itterative process.
    r=[1 1 2 2 3 3 4 4 5];
    options.maxdist=1;
    options.savememory=true;
    options.r_bins=7;
    options.a_bins=10;
    Storage=ItterativeShapeContext(Vertices1,Vertices1l,Verticesn,Verticesnl,r,options);
  
    % Show the result
    figure, 
    subplot(3,3,1), hold on; axis equal;
    plot3(Vertices1(:,2),Vertices1(:,1), Vertices1(:,3),'r.');
    plot3(Verticesn(:,2),Verticesn(:,1), Verticesn(:,3),'g.');
    for i=1:min(5,length(Storage))
        subplot(3,3,i+1), hold on;
        if(i<min(4,length(Storage))), k=i; else k=length(Storage); end
        Vertices1fit=Storage(k).Vertices;
        Matches =Storage(k).Matches;
        P1=Vertices1fit(Matches(:,1),:);
        Pn=Verticesn(Matches(:,2),:);
        FV.faces=Facesn; FV.vertices=Verticesn;
        patch(FV,'facecolor',[0 0 1],'edgecolor','none','facealpha',0.5);
        FV.faces=Faces1; FV.vertices=Vertices1fit;
        patch(FV,'facecolor',[1 0 0],'edgecolor','none','facealpha',0.5); 
        
        plot3([P1(1:5:end,1) Pn(1:5:end,1)]',[P1(1:5:end,2) Pn(1:5:end,2)]',[P1(1:5:end,3) Pn(1:5:end,3)]','g');
    end

    % Compare with manual Landmark Distance
    Cost=Storage(end).Cost;
    Vertices1fit=Storage(end).Vertices;
    Matches =Storage(end).Matches;
    clear('Storage');
    
    pp=load('images3D\land004.mat'); L1=pp.L;
    pp=load(['images3D\land' strj '.mat']); Ln=pp.L;

    % Warp the landmarks with a bspline grid constructed from the
    % corresponding points
    [~,~,L1fit]=point_registration_diff(round(size(I1)/2),Vertices1/2,Vertices1fit/2,struct('MaxRef',r(end),'Verbose',false),L1/2);
    L1fit=L1fit*2;
    
    d=sqrt(sum((L1fit-Ln).^2,2));
	disp(['landmark distance mean and std (pixels) : ' num2str(mean(d)) ' ' num2str(std(d))]);
	
    % Show the comparison
    subplot(3,3,7), hold on; axis equal;
     FV.faces=Faces1; FV.vertices=Vertices1;
     patch(FV,'facecolor',[1 0 0],'edgecolor','none','facealpha',0.5);
     plot3(L1(:,1),L1(:,2),L1(:,3),'m*','markersize',8);
     view(3);
    subplot(3,3,8),hold on; axis equal;
     FV.faces=Facesn; FV.vertices=Verticesn;
     patch(FV,'facecolor',[0 1 0],'edgecolor','none','facealpha',0.5);
     plot3(Ln(:,1),Ln(:,2),Ln(:,3),'m*','markersize',8);
     view(3);
    subplot(3,3,9), hold on; axis equal;
     FV.faces=Faces1; FV.vertices=Vertices1fit;
     patch(FV,'facecolor',[1 0 0],'edgecolor','none','facealpha',0.5);
     plot3(L1fit(:,1),L1fit(:,2),L1fit(:,3),'m*','markersize',8);
     view(3);
    drawnow
    
     % Store the contour-desription made with shape-context fitting
     Vertices=Vertices1fit; 
     Faces=Faces1; 
     Landmarks=L1fit;
     RefSurface=Verticesn;
	 
	 % Warp vertices to nearest boundary to increase Dice coefficient
	 Ic=imfill(polygon2voxel(struct('faces',Faces,'vertices',Vertices),size(In),'clamp',false),'holes');
			
	 Verticesnew=DistanceFieldCP(Vertices,Ic,[],In);
	 [~,~,Vertices]=point_registration(size(In),Vertices,Verticesnew,struct('MaxRef',r(end),'Verbose',false));

     save(['images3D\surface' strj '.mat'],'Vertices','Faces','Landmarks','RefSurface');
end


