clear all;
addpath '..\APROSFUNCTIONS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                               %%
%-----------------------------------------------------------------%
%                           File: bc3.m                           %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
%% Init source folder and target file ----------------------------%
source = 'fullmodel';            target = 'bcs.csv';              %
                                                                  %
%% Load parameter names and values from each case ----------------%
files = dir(strcat(source,'\*txt'));                              %
TIT   = {};                      TAB   = [];                      %
                                                                  %
for i = 1:length(files)                                           %
    file = files(i).name;                                         %
    file = sprintf('%s\\%s',source,file);                         %
    temp = ASignalToMatrix(file);                                 %
    % ignore timestep                                             %
    TAB(i,:) = temp(end,2:end);                                   %
                                                                  %
    if i ==1                                                      %
       TIT =  ASignalToTitle(file);                               %
       % ignore timestep                                          %
       TIT = TIT(2:end);                                          %
    end                                                           %
end                                                               %
                                                                  %
%% Write parameters for each case into one csv-file --------------%
fid = fopen(target, 'w') ;                                        %
                                                                  %
for i=1:length(TIT)                                               %
    temp = TIT(i);                                                %
    for j=1:length(temp)                                          %
        fprintf(fid, '%s,', temp{j}) ;                            %
    end                                                           %
end                                                               %
                                                                  %
fprintf(fid, '\n') ;                                              %
fclose(fid) ;                                                     %
dlmwrite(target, TAB, '-append') ;                                %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%