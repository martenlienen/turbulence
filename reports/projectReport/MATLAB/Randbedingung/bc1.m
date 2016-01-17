clear all;
addpath '..\APROSFUNCTIONS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                               %%
%-----------------------------------------------------------------%
%                           File: bc1.m                           %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
%% Create and open file to write queue to ------------------------%
AQW  = AQueueWriter('bc1-queue.txt');                             %
                                                                  %
%% Create Apros Queue Support ------------------------------------%
AQS = AQueueSupport(AQW);                                         %
                                                                  %
%% Parameter study -----------------------------------------------%
% Set point of SP47 is to be varied between 60 and 85             %
ITEM1 = AQS.getItem('SP47 SP_VALUE');                             %
steps = 60:5:85;                                                  %
                                                                  %
for i = steps                                                     %
    % Modify set point                                            %
    ITEM1.setValue(i);                                            %
    % Simulate for 40000s                                         %
    AQS.do(40000);                                                %
    % Write resuls to new initial condition                       %
    AQS.saveIC(sprintf('%d',i));                                  %
    AQW.writee();AQW.writee();AQW.writee();                       %
end                                                               %
                                                                  %
%% Close file ----------------------------------------------------%
AQW.finalize();                                                   %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%