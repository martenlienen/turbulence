function [ U,Y ] = reducedata( file, divB,divX,divY,div )
%% Function to reduce the number of supporting points
%  In ParaView you can export data along a path. You have to specify the 
%  resolution (a very high numer >1000000 to capture a cells - even
%  for streched mesh). Normally the path has more than one point in the 
%  middle of the cell. However the values (pressure, velocity...) keeps
%  constant along the cell. The data has to be reduced to one point for
%  each cell (if possible in the middle of the cell).
%
%  Input:
%    file    - file name of csv-file from paraview
%    divB    - column of data to devide by
%    divX    - column of position
%    divY    - column of data
%    div     - i.e. div=2 means: that only half of the data is considered 
%                                       (symmetry)

    % load data
    M = csvread(file,1,0);

    % data by which to divede
    U = M(1:end/div,divB);
    
    % idices of the border of each cell
    I = [];

    % find idices of the border of each cell
    temp = U(1); l = length(U); abba = 1; b=0;
    
    while(1 == abba)
        a = find(U(b+1:end)==temp,1,'first')+b;

        b = find(U(a+1:end)~=temp,1,'first');

        if isempty(b)
            b = l;
        else
            b = b + (a - 1); 
        end

        I(end+1) = a; I(end+1) = b;

        if b == l
            abba = 0;
        else
            temp = U(b+1);
        end
    end

    % reduce input data (one point for one cell)
    M = M(I,:);

    % data within cell is constant
    U = M(1:2:end,divY);
    
    % position of cell middle point
    Y = (M(1:2:end,divX)+M(2:2:end,divX))/2;

end

