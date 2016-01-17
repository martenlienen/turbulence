%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Class referancing a property of a module            %%
%-----------------------------------------------------------------%
%                       File:  AQueueItem.m                       %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
classdef AQueueItem  < handle                                     %
                                                                  %
    properties                                                    %
        name                                                      %
        AW                                                        %
    end                                                           %
                                                                  %
    methods                                                       %
                                                                  %
%% Initialize object ---------------------------------------------%
        function self = AQueueItem(AW,name)                       %
            self.name = name;                                     %
            self.AW = AW;                                         %
        end                                                       %
                                                                  %
%% Modify value of property --------------------------------------%
        function setValue(self,v)                                 %
            self.AW.writei(sprintf('MODI %s %d',self.name,v));    %
        end                                                       %
    end                                                           %
                                                                  %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%