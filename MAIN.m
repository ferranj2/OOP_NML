clc;
clear;
close all;

%Axis to visualize the fibers.
ax = custom_axis;
ax.Color = [0,0,0];
axis(ax,'equal');

%Load the fiber contours.
load('fiber_copy.mat');
fibers = length(fiber_copy);

%Define polygon objects
source_poly = polygon.CreateRegularByLength(0,0,3,1);
source_poly.SetCanvas(ax);
source_poly.Show;
source_poly.SetName('Extruder');
source_poly.SetColor([0,1,1]);

target_poly = polygon.CreateRegularByLength(0,0,3,1);
target_poly.SetCanvas(ax);
target_poly.Show;
target_poly.SetName('Offset')
target_poly.SetColor([0,1,0]);


%Cell arrays to store various nurbs.
%upper_nurbs = cell(fibers,1);
%lower_nurbs = cell(fibers,1);

offset = 0.09;
max_passes = 5;
N = 5;
lambda = 0.25;

outer = [4,7,10,11];
smallO = [5,8];
bigO = [6,9];
select = [15,20,22,28,30,35,37,40,42,43];
select2 = [14,26,29,33,36,38,41,44];%TL
select3 = [12,27,34,39];%TR

for ii = select3
    sides = length(fiber_copy{ii});
    source_poly.RedefineFromList(sides,fiber_copy{ii},true)
    if source_poly.perimeter < 0.0001
        continue;
    end%if
    source_poly.orientation
    if source_poly.orientation == 1
        mult = -1;
    else
        mult = 1;
    end%if
    %Apply Laplacian Smoothing
    for jj = 1:max_passes
        source_poly.Smooth(lambda,N);
        title(ax,['Fiber #',num2str(ii)])
    end%ii        
    stamp = Imprint(source_poly);

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

    %target_poly.CreateFromOffset(source_poly,offset);
    %target_poly.SetCanvas(ax);
    %target_poly.Show(ax);
    
    %{
    target_poly.RedefineFromList(...
        sides,...
        source_poly.Offset(offset),...
        true);
    %}
    pause(0);
end%ii