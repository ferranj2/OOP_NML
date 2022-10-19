%% Weighted B-spline distributor
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% Given a cell array of B-splines and associated weights, infers the order
% of the splines, allocates a buffer for storing polynomial coefficients of
% the superposition of the splines, and distributes the scaled piecewise
% splines over said buffer. This routine assumes that the weights provided
% correspond to the knots in the B-spline defitnitions not including the
% so-called "padding" knots. These are the potentially repeated knots
% towards the beginning or end of the defining knot sequence. Anyways, the
% weights can be the NURBS weights used to create the rational B-splines
% or the coordinates of the control points.
%% Formulae
% $R^{n}(x) = \sum_{i=1}^{k}B_{i}^{p}w_{i} = 
% \left[\begin{array}{ccccc}A_{1}&B_{1}&C_{1}&...&D_{1}\\A_{2}
% &B_{2}&C_{2}&...&D_{2}\\A_{3}&B_{3}&C_{3}&...&D_{3}\\...&...&...&...&...\\
% A_{k}&B_{k}&C_{k}&...&D_{k}\end{array}\right] \left[\begin{array}{c}1\\
% x\\x^{2}\\...\\x^{p}\end{array}\right]$
%% Required Plugins
% * none
%% Changelog
%  v1.0,(08/09/2022): Initial Release. Reinvented the wheel once again.
%% Syntax
% * INPUT(*bsp*): A cell array that contains square arrays that each
% represent a B-spline.
% * INPUT(*W*): A weight vector that will be used to scale the B-splines.
% This can be the NURBS weights (for constructing rational splines) or the
% coordinates of the control points.
% * OUTPUT(*buffer*): An array of polynomial coefficients valid between the
% implied knots that defined the B-splines.

%% Function definition
function buffer =  bspldis(bsp,W)
if iscell(bsp) == 0 %Validate list of splines
    error('B-splines must be input as a cell array.');
end
[rB,cB] = size(bsp); %Get dimensions of cell.
if rB ~=1 && cB ~= 1
    error('B-splines must be input as a column or row cell.');
end
if rB > cB
    splines = rB;
else
    splines = cB;
end
[rBold,cBold] = size(bsp{1}); %Dimension buffers now checks size of first B-spline.
if rBold ~= cBold %All B-splines form square arrays.
    error(['B-spline #',num2str(1),' does not have a square array.']);
end
for ii = 2:splines
    [rB,cB] = size(bsp{ii}); %Dimension buffers now checks size of first B-spline.
    if rBold ~= rB
        error(['Spline #',num2str(ii),'has different rows from spline #',num2str(ii-1)]);
    end
    if cBold ~= cB
        error(['Spline #',num2str(ii),'has different columns from spline #',num2str(ii-1)]);
    end
    rBold = rB;
    cBold = cB;
end
p = rB - 1; %Infer the order of the B-splines from inspecting the arrays.

%Validate weight input.
[rw,cw] = size(W);
if rw ~=1 && cw ~= 1
    error('Weights must be supplied as a column or row vector.');
end
if rw > cw
    knots = rw;
else
    knots = cw;
end
%Allocate output buffer;
buffer = zeros(knots-1,p+1); %Allocate space for the polynomial coefficients.

%Begin scaled spline distribution.
for ii = 1:p %Distribute first "p" splines.
    bsp{ii}
    %for jj = (p - ii + 2):(p + 1) %Splines of correct order.
    for jj = 1:ii %Splines of correct order. 1

        %jj
        %p + 1 + jj-ii
        
        %p + 2 - ii: p + 1 
        for ww = 1:(p + 1) %Polynomial coefficients.
            %buffer(jj - p + ii - 1,ww) = buffer(jj - p + ii - 1,ww) + bsp{ii}(jj,ww)*W(1); %Don't forget to scale.
            buffer(jj,ww) = buffer(jj,ww) + bsp{ii}(p + 1 + jj - ii,ww)*W(1); %Don't forget to scale.

        end
    end
end
for ii = (p + 1):(splines - p) %All splines/intervals that do not involve a "ghost" knot.
    for jj = (ii - p):ii %ii is only offset by "p" over this spline range.
        for ww = 1:(p + 1) %Polynomial coefficients.
            buffer(jj,ww) = buffer(jj,ww) + bsp{ii}(jj - ii + p + 1,ww)*W(ii); %Don't forget to scale.
        end
    end
end
segments = knots - 1;
for ii = (splines - p + 1):splines %Right-end "ghosts."
    for jj = (segments -(splines - ii)):(segments) %Index of the registers.
        for ww = 1:(p + 1) %Polynomial coefficients.
           buffer(jj,ww) = buffer(jj,ww) + bsp{ii}(jj + splines - ii - segments + 1,ww)*W(knots); %Don't forget to scale.
        end
    end
end
end