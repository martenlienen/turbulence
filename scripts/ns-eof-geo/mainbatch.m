global obstacle file folder

%% location of the configuaration-file
folder = '.';
file   = 'conf_channel_praesi.xml';

%% type of obstacle
obstacle = obstaclecircle(0.2,0.48,0.35);
% obstacle = obstaclesquare(0.2,0.49,0.15);
% obstacle = obstacleplate(0.2,0.49,0.40,0.05,25);

%% perfom task (create .geo file)
main();
