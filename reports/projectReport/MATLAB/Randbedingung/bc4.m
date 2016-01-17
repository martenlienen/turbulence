clear all;
addpath '..\APROSFUNCTIONS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                               %%
%-----------------------------------------------------------------%
%                           File: bc4.m                           %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
% create file to write queue to                                   %
AQW  = AQueueWriter('bc4-queue.txt');                             %
AQS = AQueueSupport(AQW);                                         %
                                                                  %
ITEMS = AQS.getItems('bcs.csv');                                  %
ITEMS.renameProperty('SP47 SP_VALUE','XA01 ANALOG_VALUE');        %
ITEMS.setValue(70);                                               %
                                                                  %
%% Close file ----------------------------------------------------%
AQW.finalize();                                                   %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%