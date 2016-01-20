clear all
clc
set(0,'RecursionLimit',100000)

%% number of cells in x and y direction
numberOfPointsX = 128;
numberOfPointsY = 128;

%% setup Domain

figure(2)
hold on

%img = imread('airfoil.png');             %# Load a sample image
%h = image([0 numberOfPointsX],[numberOfPointsY 0],img);        %# Plot the image

[topoGrid, coordsX, coordsY] = getDomain(numberOfPointsX,numberOfPointsY);
seedPoint = round(ginput(1));
topoGrid = FloodFill(topoGrid, seedPoint);




%% Write obstacles to file

% specifiy filename
fileID = fopen('small-geo.txt','w');

for i=1:numberOfPointsX
    for j=1:numberOfPointsY
        if(topoGrid(i,j)==1 || topoGrid(i,j)==2)
           fprintf(fileID,'%5d %5d %5d\n',i-1,j-1,0); 
        end
    end
end

%h.Visible = 'off';
