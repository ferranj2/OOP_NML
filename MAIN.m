clc;
clear;
close all;


%Axis to visualize the fibers.
ax = custom_axis;
ax.Color = [0,0,0];
ax.DataAspectRatioMode = 'manual';
axis(ax,'equal');

%Load the fiber contours.
load('fiber_copy.mat');
fibers = length(fiber_copy);
path_nurbs = cell(fibers,1); %This will store NURBS paths.
upper_nurbs = cell(fibers,1); %These are offset from the path NURBS.
lower_nurbs = cell(fibers,1); %These are offset from the
p = 2;
g = 3;

%Define an "Extruder" polygon.
source_poly = polygon.CreateRegularByLength(0,0,3,1);
source_poly.SetRefreshRate(4);
source_poly.SetCanvas(ax);
source_poly.Show;
source_poly.ToggleAABB;
source_poly.SetName('Extruder');
source_poly.SetColor([0,1,1]);

%Define an "Offset Polygon.
target_poly = polygon.CreateRegularByLength(0,0,3,1);
target_poly.SetRefreshRate(4);
target_poly.SetCanvas(ax);
target_poly.Show;
target_poly.ToggleAABB;
target_poly.SetName('Offset')
target_poly.SetColor([0,1,0]);

%Define a "Path" polygon
path_poly = polygon.CreateRegularByLength(0,0,3,1);
path_poly.SetCanvas(ax);
path_poly.Show;
path_poly.SetName('Path');
path_poly.SetColor([1,1,1]);

%Cell arrays to store various nurbs.


offset = 0.09;
max_passes = 2;
N = 5;
lambda = 0.25;
Rf = 0.04; %fillet_radius 

outer = [4,7,10,11];
smallO = [5,8];
bigO = [6,9];
select = [15,20,22,28,30,35,37,40,42,43];
select2 = [14,16,21,26,29,33,36,38,41,44];%TL
select3 = [12,23,27,32,34,39];%TR
select4 = [24,25,31];
all = 1:fibers;



%select4 has some issues
%25[mm]
%Add xy labels with units.
[select,select2,select3,select4]
for ii = smallO
    %for ii = select1(2:end)
    sides = length(fiber_copy{ii});
    source_poly.RedefineFromList(sides,fiber_copy{ii},true)
    if ii == select(1)
        %input('Ready?')
    end%if
    if source_poly.perimeter < 0.0001
        continue;
    end%if
    
    if source_poly.open
        if source_poly.orientation == 1
            mult = -1; % + 1
        else
            mult = +1; %-1
        end%if
    else
        mult = 1;
    end
    
    
    
    
    %Apply Laplacian Smoothing
    for jj = 1:max_passes
        source_poly.Smooth(lambda,N);
        title(ax,['Fiber #',num2str(ii)])
    end%ii  
    %source_poly.aabb.FrameCanvas;

    source_poly.DetectSharpTurns(Rf);
    %stamp = Imprint(source_poly);

    %The "target" polygon data structure is redefined as the offset of the
    %"source" polygon data structure. 
    target_poly.RedefineFromList(...
        source_poly.sides,...
        source_poly.Offset(mult*offset),...
        source_poly.open);
    if ~target_poly.simple
        %if self-intersections detected.
    end

    stamp2 = Imprint(target_poly);

    %The "path" polygon will define the centerline nurbs.
    path_poly.RedefineAsWeightedAverage(source_poly,1,target_poly,1);
    
    %Fit NURBS through the path polygon.
    path_nurbs{ii} = nurbs.CreateFromPolygon(path_poly,p,g);
    path_nurbs{ii}.SetCanvas(ax);
    path_nurbs{ii}.Show;
    path_nurbs{ii}.sketches.Curve.LineStyle = '--';
    
    %Offset in the "+" direction.
    upper_nurbs{ii} = nurbs.CreateFromNURBSOffset(path_nurbs{ii},offset*0.5);
    upper_nurbs{ii}.SetCanvas(ax);
    upper_nurbs{ii}.SetColor([0,0,1]); 
    upper_nurbs{ii}.Show;
    
    %Offset in the "-" direction.
    lower_nurbs{ii} = nurbs.CreateFromNURBSOffset(path_nurbs{ii},-offset*0.5);
    lower_nurbs{ii}.SetCanvas(ax);
    lower_nurbs{ii}.SetColor([1,0,0]); 
    lower_nurbs{ii}.Show;
        
    pause(0);
end%ii
