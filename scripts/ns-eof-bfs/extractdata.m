clear all;

%% Skript to load from VTK-files the values along a path for every time step

%% Setup

% x-positon [m]
xpositionm = 0.35;

% folder with vtk-files
folder = 'karmanfree';

% name of the csv-file, to which to write data
resfile = 'karmanfree350.csv';

% numbers of domains to consider (in right order)
dekompy = [0];

% start string of velocity section
vsection = 'VECTORS velocity float';

%% List all vtk-files and sort them after iteration number -> index

% list all files in folder
files = dir(folder);

% delete '.' and '..'
files = files(3:end);

% number of files per timestep 
fnr = length(files)/length(dekompy);

index = [];

for i = 1:fnr 
    % read iteration numer
   temp = strsplit(files(i).name,'.');
   index(end+1,1) = i; 
   index(end  ,2) = str2double(temp(end-1));
end

% sort
index = sortrows(index,2);



%% Load velocity v for each cell at every timestep
%  MMMM is created (row: timesteps, cols: position in y-direction)

offsetj = 2;

tic
for dekomp = dekompy
sprintf('outer-loop = %d', dekomp )
    for i = 1:fnr
        sprintf('inner-loop = %d', i )
        filename = sprintf('%s\\%s',folder,char(files(index(i,1)).name));

        if(strfind(filename,'.0.0.'))
            filename = strrep(filename,'.0.0.',sprintf('.%d.0.',dekomp));
        else
            filename = strrep(filename,'.0.',sprintf('.%d.',dekomp));
        end

        if i == 1
            % first line of velocity section
            vp = findline(filename,vsection);

            % numerb of nodes in each direction
            ttp =findline(filename,'DATASET STRUCTURED_GRID') + 1;
            tp = getline(filename, ttp);
            dimensions = sscanf(tp,'DIMENSIONS %d %d %d');
            
            %% init MMM
            if offsetj==2
                MMMM = zeros(1+fnr, length(dekompy)*(dimensions(2)-1)+1);
            end
            
            % transform x-position from [m] in i
            M = dlmread(filename, ' ',[ttp+1 0 ttp+dimensions(1) 0]);
            xposition = find(M>=xpositionm,1,'first');

            % get y-coordinates
            M = dlmread(filename, ' ',[ttp+1 0 ttp+dimensions(1)*dimensions(2) 2]);
            M = M(xposition:dimensions(1):end,2);

            M = (M(1:end-1)+M(2:end))/2;

            % write y-coordinates
            MMMM(1,offsetj : offsetj + dimensions(2) - 2) = M';
        end

        % load velocities
        M = dlmread(filename, ' ',[vp 0 vp+(dimensions(1)-1)*(dimensions(2)-1)-1 2]);
        % delete v and w & take values along path
        M = M(xposition:(dimensions(1)-1):end,1);

        if i == 1 && offsetj == 2
            % write timestep
            MMMM(i + 1, 1) = index(i,2);
        end
        
        % write v-values
        MMMM(i + 1, offsetj : offsetj + dimensions(2) - 2) = M';
    end

    offsetj = offsetj + dimensions(2) - 1;
end
toc


% dummy-plot
for i = 2:size(MMMM,1)
   plot(MMMM(i,2:end),MMMM(1,2:end));hold on;grid on; 
end

% write to CSV-File
csvwrite(resfile,MMMM);
