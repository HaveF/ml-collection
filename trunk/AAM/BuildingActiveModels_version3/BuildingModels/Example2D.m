% This is an Example of using Shape-Context to create a corresponding
% point model of a set of 10 photos (2D) of hands.
%
% It loads 10 datasets and convert them to a point
% description of the outer contour. ShapeContext is used to
% iteratively register one of the datasets to all other datasets.
% The warping between datasets is kept diffeomorphic, by constraining
% the jacobian of the b-spline transformation grid.
%
% To use the constructed point model to create an ASM or AAM
% run the examples in ActiveModels_version7
%
% Example is written by D.Kroon University of Twente 
% March 2011 updated January 2012


% Add all the needed folders and functions to Matlab search path
functionname='Example2D.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'functions'])
addpath([functiondir 'isocontour_version2'])
addpath([functiondir 'ShapeContext_version2'])
addpath([functiondir 'nonrigid_version25'])
% This will compile the c-code to fast multi-threaded mex-files. The
% example also works without the compile c-files but will be slower.
cd([functiondir 'ShapeContext_version2']);
try
    compile_c_files,
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
k=4;
strk=['000' num2str(k)]; strk=strk(end-2:end);
disp(['Master : ' strk]); drawnow;
filename=['images2D\segm' strk '.png'];

% Convert the segmented image to contour data
[I1,Lines1,Vertices1,Lines1l,Vertices1l]=GetImageContour(filename);

% Show the dataset
figure, imshow(I1), hold on;
V1=Vertices1(Lines1(:,1),:);
V2=Vertices1(Lines1(:,2),:);
plot([V1(:,2) V2(:,2)]',[V1(:,1) V2(:,1)]','b');
plot(Vertices1(:,2),Vertices1(:,1),'r.');
drawnow

for j=1:10,
    % Load a dataset
    strj=['000' num2str(j)]; strj=strj(end-2:end);
    filename=['images2D\segm' strj '.png'];
    disp(['Load Dataset : ' filename]); drawnow('expose');
    % Convert the segmented image to contour data
    [In,Linesn,Verticesn,Linesnl,Verticesnl]=GetImageContour(filename);
    
    % Itterative register the contours using shape context
    % regulize the matching with a b-spline grid which is refined from
    % coarse to fine during this itterative process.
    r=[1 1 1 2 2 2 3 3 3 4 4 5];
    options.maxdist=1;
    options.savememory=true;
    Storage=ItterativeShapeContext(Vertices1,Vertices1l,Verticesn,Verticesnl,r,options,Lines1,Linesn);
    
    % Show the result
    figure,
    subplot(3,3,1), hold on; axis equal;
    plot(Vertices1(:,2),Vertices1(:,1),'r.');
    plot(Verticesn(:,2),Verticesn(:,1),'g.');
    for i=1:5
        subplot(3,3,i+1), hold on;
        if(i<4), k=i; else k=length(Storage); end
        Cost=Storage(k).Cost;
        Vertices1fit=Storage(k).Vertices;
        Matches =Storage(k).Matches;
        P1=Vertices1fit(Matches(:,1),:);
        Pn=Verticesn(Matches(:,2),:);
        plot(Vertices1fit(:,2),Vertices1fit(:,1),'b.');
        plot(Verticesn(:,2),Verticesn(:,1),'r.');
        Cost=Cost./max(max(Cost(:),eps));
        for t=1:size(Matches,1)
            C = Cost(t)*5;
            plot([P1(t,2) Pn(t,2)],[P1(t,1) Pn(t,1)],'g','LineWidth',C+0.2);
        end
    end
    
    % Compare with manual Landmark Distance
    Cost=Storage(end).Cost;
    Vertices1fit=Storage(end).Vertices;
    Matches =Storage(end).Matches;
    
    pp=load(['images2D\land' strk '.mat']); L1=pp.L;
    pp=load(['images2D\land' strj '.mat']); Ln=pp.L;
    
    % Warp the landmarks with a bspline grid constructed from the
    % corresponding points
    [~,~,L1fit]=point_registration_diff(size(I1),Vertices1,Vertices1fit,struct('MaxRef',r(end),'Verbose',false),L1);

    % Show the comparison
    subplot(3,3,7), hold on;
    plot(Vertices1(:,2),Vertices1(:,1),'r.');
    plot(L1(:,2),L1(:,1),'m*','markersize',8);
    subplot(3,3,8),hold on;
    plot(Verticesn(:,2),Verticesn(:,1),'g.');
    plot(Ln(:,2),Ln(:,1),'c*','markersize',8);
    subplot(3,3,9), hold on;
    plot(Verticesn(:,2),Verticesn(:,1),'r.');
    plot(Ln(:,2),Ln(:,1),'c*','markersize',8);
    plot(Vertices1fit(:,2),Vertices1fit(:,1),'b.');
    plot(L1fit(:,2),L1fit(:,1),'m*','markersize',8);
    drawnow
    
    % Store the contour-desription made with shape-context fitting
    Vertices=Vertices1fit;
    Verticesl=Storage(i).Verticesl;
    Lines=Lines1;
    Linesl=Lines1l;
    Landmarks=L1fit;
    RefContour=Verticesn;
    save(['images2D\contour' strj '.mat'],'Vertices','Lines','Landmarks','RefContour');
end

% Show all results
figure
for j=1:10,
    % Load a dataset
    strj=['000' num2str(j)]; strj=strj(end-2:end);
    filename=['images2D\train' strj '.jpg'];
    load(['images2D\contour' strj '.mat']);
    I=imread(filename);
    subplot(3,4,j),imshow(I);  hold on;
    V1=Vertices(Lines(:,1),:); V2=Vertices(Lines(:,2),:);
    plot(RefContour(:,2),RefContour(:,1),'g.');
    plot([V1(:,2) V2(:,2)]',[V1(:,1) V2(:,1)]','b');
    plot(Landmarks(:,2),Landmarks(:,1),'m*','markersize',4);
end

