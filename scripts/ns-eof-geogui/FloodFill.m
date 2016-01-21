function [ topoGrid ] = FloodFill( topoGrid, seedPoint )
%FLOODFILL Summary of this function goes here
%   Detailed explanation goes here

%% recursive algorithm to detect all obstacle cells

%% setup variables
boundary = 2;
unset = -5;
set = 1;
limitX = size(topoGrid,1);
limitY = size(topoGrid,2);

%disp(seedPoint);


%% check if seedpoint is in boundary
while topoGrid(seedPoint(1),seedPoint(2))==boundary
    warning('Please choose another seedpoint')
    seedPoint = round(ginput(1));
end

%% check all directions
if topoGrid(seedPoint(1),seedPoint(2))==unset
    topoGrid(seedPoint(1),seedPoint(2))=set;
end

if topoGrid(seedPoint(1),seedPoint(2)+1)<=limitY && topoGrid(seedPoint(1),seedPoint(2)+1)==unset
    topoGrid = FloodFill(topoGrid,[seedPoint(1) seedPoint(2)+1]);
end

if topoGrid(seedPoint(1),seedPoint(2)-1)<=limitY && topoGrid(seedPoint(1),seedPoint(2)-1)==unset
    topoGrid = FloodFill(topoGrid,[seedPoint(1) seedPoint(2)-1]);
end

if topoGrid(seedPoint(1)+1,seedPoint(2))<=limitX && topoGrid(seedPoint(1)+1,seedPoint(2))==unset
    topoGrid = FloodFill(topoGrid,[seedPoint(1)+1 seedPoint(2)]);
end

if topoGrid(seedPoint(1)-1,seedPoint(2))<=limitX && topoGrid(seedPoint(1)-1,seedPoint(2))==unset
    topoGrid = FloodFill(topoGrid,[seedPoint(1)-1 seedPoint(2)]);
end

