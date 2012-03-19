function Storage=ItterativeShapeContext(Vertices1,Vertices1l,Verticesn,Verticesnl,r,options,Lines1,Linesn)
%
% Storage = ItterativeShapeContext(P1,P2,PL1,PL2,S,R,Options)
%
%  inputs,
%    P1 : Point List Nx2 or Nx3 describing the object contour in
%            the first dataset
%    P1 : Point List Mx2 or Mx3 describing the object contour in
%            the second dataset
%    PL1 : Set equal to P1, or to a point list with more samples describing
%         the shame shape as P1
%    PL2 : Set equal to P2, or to a point list with more samples describing
%         the shame shape as P2
%	 R : A vector with the number of transformattion grid refinements for each step
%
%    Options : A struct with options:
%      Options.r_bins : Number of Log Distance bins in the
%                           matching histogram(default 15)
%      Options.a_bins : Number of Angle bins in the
%                           matching histogram (default 15)
%      Options.rotate : Rotation Normalization of the points using
%                       eigen vector analyis 0,1 or 2.
%                       0: No correction (default), 1: Minimal align
%                       rotation(1), 2: Rotation with Heavyside flip.
%      Options.method : The matching method used:
%                       0 : Munkres One to one matching
%                       1 : Minimum One to multiple matching (default)
%      Options.maxdist : The maximum matching distance between normalized
%                        points (default 5). The Point-sets are normalized
%                        to have an average radius of one to their mean value.
%
%  outputs,
%    Storage : A struct with all the matchting information such as the b-spline
%				transformation grid
%
% See also EXAMPLE2D, EXAMPLE3D, SHAPECONTEXT
%
% This function is written by D.Kroon University of Twente (March 2011)
Vertices1fit=Vertices1;
Vertices1lfit=Vertices1l;
Storage=struct;
if(nargin>7)
    texturesize=ceil(max([max(Verticesn,[],1);max(Vertices1,[],1)],[],1));
    In=drawObject(Verticesn,texturesize,Linesn)>0;
end

options2=struct;
if(isfield(options,'maxdist')),options2.maxdist=options.maxdist; end
if(isfield(options,'r_bins')), options2.r_bins=options.r_bins; end
if(isfield(options,'a_bins')), options2.a_bins=options.a_bins; end

options3=struct;
if(isfield(options,'savememory')), savememory=options.savememory; else savememory=false; end
    
for i=1:length(r)
    options2.rotate=i<3;
    % Detect matching vertices, using shape context as feature
    
    [Matches,Cost]=ShapeContext(Vertices1fit,Verticesn,Vertices1lfit,Verticesnl,options2);
    
    % Fit a bspline grid to transform points Vertices1 into Verticesn
    options3.Verbose=false;
    options3.MaxRef=r(i);

    % Calculate the b-spline warped-points positions
    [Vertices1fit,Vertices1lfit]=WarpVertices(Vertices1,Verticesn,Vertices1l,Matches,options3,savememory);

    % Calculate distance between matches
    E=mean(sqrt(sum((Vertices1fit(Matches(:,1),:)-Verticesn(Matches(:,2),:)).^2,2)));
    if(nargin>7)
        I1=drawObject(Vertices1fit,texturesize,Lines1)>0;
        dice=2*sum(In(:)&I1(:))/(sum(In(:))+sum(I1(:)));
        disp(['Shape Context Itteration: ' num2str(i) ' Refinements Grid : ' num2str(r(i)) ' Mean Point Distance : ' num2str(E) ' Dice Coeff : ' num2str(dice)]);
    else
        disp(['Shape Context Itteration: ' num2str(i) ' Refinements Grid : ' num2str(r(i)) ' Mean Point Distance : ' num2str(E)]);
    end
    
    % Store information about the previous transformation field and
    % shape context matching
    Storage(i).Vertices=Vertices1fit;
    Storage(i).Verticesl=Vertices1lfit;
    Storage(i).Matches=Matches;
    Storage(i).Cost=Cost;
    Storage(i).E=E;
end


function [Vertices1fit,Vertices1lfit]=WarpVertices(Vertices1,Verticesn,Vertices1l,Matches,options,savememory)
if(savememory), scale=2; else scale=1; end
% 
offs=5;

% Scale down to save memory
Vertices1=Vertices1/scale;
Vertices1l=Vertices1l/scale;
Verticesn=Verticesn/scale;

% Boundaries of vertices + padding
minb=min(min(Vertices1,[],1),min(Verticesn,[],1))-offs;

% Substract boundary
Vertices1=bsxfun(@minus,Vertices1,minb);
Verticesn=bsxfun(@minus,Verticesn,minb);
Vertices1l=bsxfun(@minus,Vertices1l,minb);

% size of volume
Isize=ceil(max(max(Vertices1,[],1),max(Verticesn,[],1)))+offs;

% The Matched points
P1=Vertices1(Matches(:,1),:);
Pn=Verticesn(Matches(:,2),:);

[~,~,Vertices1fit,Vertices1lfit]=point_registration_diff(Isize,P1,Pn,options,Vertices1,Vertices1l);

% Re-add boundary
Vertices1fit=bsxfun(@plus,Vertices1fit,minb);
Vertices1lfit=bsxfun(@plus,Vertices1lfit,minb);

% Scale up to save memory
Vertices1fit=Vertices1fit*scale;
Vertices1lfit=Vertices1lfit*scale;



