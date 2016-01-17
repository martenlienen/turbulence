%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                   Class creating the queue                    %%
%-----------------------------------------------------------------%
%                      File: AQueueSupport.m                      %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
classdef AQueueSupport  < handle                                  %
                                                                  %
    properties                                                    %
        AW                                                        %
    end                                                           %
                                                                  %
    methods                                                       %
                                                                  %
%% Initialize object ---------------------------------------------%
        function self = AQueueSupport(AW)                         %
            % AW: a writer (i.e. AQueueWriter)                    %
            self.AW = AW;                                         %
        end                                                       %
                                                                  %
%% Returns reference to a modifiable property --------------------%
        function IT = getItem(self,item)                          %
            IT = AQueueItem(self.AW,item);                        %
        end                                                       %
                                                                  %
%% Returns reference to a modifiable property --------------------%
        function IOSS = getIOSetSupport(self)                     %
            IOSS = AIOSetSupport(self,self.AW);                   %
        end                                                       %
                                                                  %
%% Save IC -------------------------------------------------------%
        function saveIC(self,ic)                                  %
            self.AW.writei(sprintf('save_ic %s',ic));             %
        end                                                       %
                                                                  % 
%% Load IC -------------------------------------------------------%
        function loadIC(self,ic)                                  %
            self.AW.writei(sprintf('load_ic %s',ic));             %
        end                                                       %
                                                                  %
%% Simulate n seconds --------------------------------------------%
        function do(self,n)                                       %
            self.AW.writei(sprintf('do %d',n));                   %
        end                                                       %
                                                                  %
    end                                                           %
                                                                  %  
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
