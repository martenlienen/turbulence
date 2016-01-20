%% Load packages
addpath('xmlfunctions');

%% Load configuration-file
conf = loadconfig(char(folder),char(file));

%% Create mesh

[x,y] = createmesh(conf);

%% Create obstacles   
z = ones(length(x),length(y));
z = obstacle(conf,x,y,z);

%% Visualize mesh and obstacles
figure 
subplot(2,1,1)
h = pcolor(x,y,z'); 
colormap(gray(2))
axis equal tight

subplot(2,1,2)
h = pcolor(x,y,z'); 
colormap(gray(2))
axis equal tight
set(h, 'EdgeColor', 'none');

set(gcf,'color','w');

%% Write obstacles to file
fileID = fopen(sprintf('%s\\%s.geo',folder,strrep(file,'.xml','')),'w');

for i=1:length(x)
    for j=1:length(y)
        if(z(i,j)==0)
           fprintf(fileID,'%5d %5d %5d\n',i-1,j-1,0); 
        end
    end
end

fclose(fileID);
