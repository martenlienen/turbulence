function [ topoGrid, coordsX, coordsY ] = getDomain( numberOfPointsX, numberOfPointsY )
%GETDOMAIN Summary of this function goes here
%   Detailed explanation goes here

%% setup variables
unset = -5;

%% create topoGrid
topoGrid = ones(numberOfPointsX, numberOfPointsY).*unset;

%% create coords-Matrices
[coordsX, coordsY] = meshgrid(1:1:numberOfPointsX, 1:1:numberOfPointsY);
coordsX = coordsX';
coordsY = coordsY';

%% prepare figure

axis([1, numberOfPointsX, 1, numberOfPointsY])
grid on
point1 = round(ginput(1));
initialPoint = point1;

%% detect user specified points and calculate line between them
while 1
    point2 = round(ginput(1));
    [topoGrid, coordsX, coordsY] = bresenham(point1, point2, topoGrid, coordsX, coordsY);
    if point2 == initialPoint
        break
    end
    point1 = point2;
end


end

