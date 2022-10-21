clc;
clear;
close all;


%Axis to visualize the fibers.
ax = custom_axis;
ax.Color = [0,0,0];
ax.DataAspectRatioMode = 'manual';
axis(ax,'equal');

%Load the fiber contours.
%load('fiber_copy.mat');
load('LVL10.mat');
fibers = length(fiber_copy);
path_nurbs = cell(fibers,1); %This will store NURBS paths.
upper_nurbs = cell(fibers,1); %These are offset from the path NURBS.
lower_nurbs = cell(fibers,1); %These are offset from the
p = 2;
g = 2;

holes = 3;
HOLES = cell(holes,1);
HOLES{1} = circle.Create3P(...
    14.343007612499999,3.205434852010000,...
    14.112008032200000,2.616045693290000,...
    13.677080966500000,2.763968715690000);
HOLES{2} = circle.Create3P(...
    5.601220977360000,13.030954240000000,...
    6.367368107950000,13.158044077400000,...
    5.794524416330000,13.343175686900000);
HOLES{3} = circle.Create3P(...
    1.710843897360000,13.276321824900000,...
    2.399953947350000,13.002428043800000,...
    1.620686252150000,13.126827021600000);
for ii = 1:holes
    HOLES{ii}.SetCanvas(ax);
    HOLES{ii}.SetColor([1,1,1]);
    HOLES{ii}.sketches.Curve.LineWidth = 1;
    HOLES{ii}.Show;
end%ii

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


offset = 0.1;
max_passes = 2; %Number of Laplacian smooths 
N = 5; %Number of Stencil points to use in the Laplacian Smooths.
lambda = 0.25; %Magnitude of the smooths
Rf = offset/2; %fillet_radius for rollercoaster.

outer = [4,7,10,11];
smallO = [5,8];
bigO = [6,9];
select = [15,20,22,28,30,35,37,40,42,43];
select2 = [14,16,21,26,29,33,36,38,41,44];%TL
select3 = [12,23,27,32,34,39];%TR
select4 = [24,25,31];
all = 1:fibers;


front1 = [1,4,7]; %L-shape
front2 = [2,5,8]; %Small O
front3 = [3,6,9]; %Big O
front4 = [10,12]; %L-shape split by holes towards bottom.
front5 = [11]; %L-shape split by holes towards top.
front6 = [13,22,38,57]; %Top-right "tear"
front7 = [14]; %Bottom-left triangle of front 5.
front8 = [15,24]; %Top-Left triangle of front 5.
front9 = [17,26]; %Top-left "tear".
front10 = [18,45,58,64,69]; %Bottom-right hole shroud (front 4).
front11 = [20,52,59,65,70]; %Top-left Triangle of front 4.
front12 = [21,53]; %Bottom-right triangle of front 4.
bad = [25,26:37,39:44,46:51,54:56,60:63];

%select4 has some issues
%25[mm]
%Add xy labels with units.
%[select,select2,select3,select4]
for ii = 34:fibers
    %for ii = select1(2:end)
    sides = length(fiber_copy{ii});
    if sides < 6
        continue;
    end%if
    source_poly.RedefineFromList(sides,fiber_copy{ii},true)
        
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
    if path_poly.perimeter < 2.5
        continue;
    end%if
    
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
       
    %input('Next?');
    pause(0);
end%ii

%{
L = zeros(fibers,1);
for ii = 1:fibers
    L(ii) = length(fiber_copy{ii});
end
%}


%{
Toggle NURBS on/off
for ii = 1:44
    upper_nurbs{ii}.sketches.Curve.Visible = 'on';
    lower_nurbs{ii}.sketches.Curve.Visible = 'on';
end
%}
control_polygons = cell(fibers,1);
for ii = 1:fibers
    if ~isempty(path_nurbs{ii})
        control_polygons{ii} = path_nurbs{ii}.CP;
    end%if
end%ii
save('control_polygons.mat','control_polygons')

