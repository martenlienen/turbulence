clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                               %%
%-----------------------------------------------------------------%
%                           File: bc2.m                           %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
%% Create and open file to write queue to ------------------------%
AQW  = AQueueWriter('bc2-queue.txt');                             %
                                                                  %
%% Create Apros Queue Support ------------------------------------%
AQS = AQueueSupport(AQW);                                         %
                                                                  %
%% Parameters of interest ----------------------------------------%
ALL = {'PO11_PRESSURE','PO11_TEMPERATURE','PO11_MIX_ENTH',...     %
    'PO11_LIQ_ENTH','PO11_STEAM_ENTH','PO11_VOID',...             %
    'PO11_NONC_MASS_FRAC'};                                       %
                                                                  %
%% Parameter study -----------------------------------------------%
% Set point of SP47 is to be varied between 60 and 85             %
steps = 60:5:85;                                                  %
                                                                  %
for i = steps                                                     %
    AQS.loadIC(sprintf('%d',i));                                  %
                                                                  %
    IOS = AQS.getIOSetSupport(); IOS.init('IO04');                %
    IOS.timestep(0.05);                                           %
                                                                  %
    IOS.monitor('SP47',{'SP_VALUE'});                             %
    IOS.monitor('PO255',ALL); IOS.monitor('PO202',ALL);           %
    IOS.monitor('PO157',ALL); IOS.monitor('PO29' ,ALL);           %
    IOS.monitor('PO204',ALL); IOS.monitor('PO267',ALL);           %
    IOS.monitor('PO34' ,ALL); IOS.monitor('PO214',ALL);           %
    IOS.monitor('PO236',ALL); IOS.monitor('PO153',ALL);           %
    IOS.monitor('PO158',ALL); IOS.monitor('PO301',ALL);           %
    IOS.monitor('PO279',ALL);                                     %
                                                                  %
    IOS.title(sprintf('res%f.txt',i));                            %
    IOS.runAndWrite(1);                                           %
                                                                  %
    AQW.writee();AQW.writee();AQW.writee();                       %
    IOS.finalize();                                               %
    AQW.writee();AQW.writee();AQW.writee();                       %
end                                                               %
                                                                  %
%% Close file ----------------------------------------------------%
AQW.finalize();                                                   %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

