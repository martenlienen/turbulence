%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                               %%
%-----------------------------------------------------------------%
%                           File: gf1.m                           %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
%% Initialize ----------------------------------------------------%
                                                                  %
% Folder to read files from                                       %
folder = 'res_stau_85_fm6';                                       %
                                                                  %
%% Process data from files ---------------------------------------%
files = dir(sprintf('%s\\*txt',folder));  RES   = [];             %
                                                                  %
for i = 1:length(files)                                           %
    file = files(i).name;                                         %
    file = sprintf('%s\\%s',folder,file);                         %
    [f0,fin,fout,w,a,p,c,e1,e2] = ...                             %
          SignalToBodeData(file,@NameToTimeperiode,2,3);          %
                                                                  %
    RES(i,:)=[f0,fin,fout,w,a,p,c,e1,e2];                         %
end                                                               %
                                                                  %
% Sort rows after the frequency                                   %
RES = sortrows(RES,1);                                            %
                                                                  %
%% Visualize data from files -------------------------------------%
BodeDataToBodeDiagram(RES(:,1),RES(:,5),RES(:,6))                 %
BodeDataToNyquistDiagram(RES(:,1),RES(:,5),RES(:,6))              %
BodeDataToError(RES(:,1),RES(:,8),RES(:,9))                       %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%