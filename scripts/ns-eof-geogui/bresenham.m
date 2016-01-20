function [ topoGrid, coordsX, coordsY ] = bresenham( point1, point2, topoGrid, coordsX, coordsY )
%BRESENHAM Summary of this function goes here
%   Detailed explanation goes here


%%% bresenham algorithm to determine a line between two points on a raster

dist = point2-point1;
direction = sign(dist);

%% determine step (direction) in which "quadrant direction" the vector
%  between point1 and point2 points
if abs(dist(1))>=abs(dist(2))
    longStep = [direction(1) 0];
    shortStep = [0 direction(2)];
    m = dist(2)/dist(1);
    t = dist(1);
else
    longStep = [0 direction(2)];
    shortStep = [direction(1) 0];
    m = dist(1)/dist(2);
    t = dist(2);
end

%% start from point 1 and walk single steps up to point 2
cPoint = point1;
error = 0;
figure(2)
hold on
for i = 1:1:abs(t)
    cPoint = cPoint + longStep;
    error = error + m;
    if abs(error)>=0.5
        cPoint = cPoint + shortStep;
        error = error - sign(m);
    end
    topoGrid(cPoint(1), cPoint(2)) = 2;
    coordsX(cPoint(1), cPoint(2)) = coordsX(cPoint(1), cPoint(2)) + error*longStep(2);
    coordsY(cPoint(1), cPoint(2)) = coordsY(cPoint(1), cPoint(2)) + error*longStep(1);
    plot(coordsX(cPoint(1), cPoint(2)), coordsY(cPoint(1), cPoint(2)),'o-');
end

