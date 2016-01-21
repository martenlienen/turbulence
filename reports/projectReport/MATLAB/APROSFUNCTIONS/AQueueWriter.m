%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%              Class responible for the Queue-file              %%
%-----------------------------------------------------------------%
%                      File:  AQueueWriter.m                      %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
classdef AQueueWriter                                             %
                                                                  %
    properties                                                    %
        FILE                                                      %
    end                                                           %
                                                                  %
    methods                                                       %
                                                                  %
%% Initialize object ---------------------------------------------%
        function self = AQueueWriter(FILE)                        %
            % Create and open file                                %
            self.FILE = fopen(FILE,'w');                          %
        end                                                       %
                                                                  %
%% Write line to file --------------------------------------------%
        function writei(self,text)                                %
            fprintf( self.FILE,'%s\n',text);                      %
        end                                                       %
                                                                  %
%% Write lines to file -------------------------------------------%
        function writeis(self,text)                               %
            l = length(text);                                     %
            if l > 1                                              %
                for t = 1:length(text)-1                          %
                    fprintf(self.FILE,'%s - \n',char(text(t)));   %
                end                                               %
            end                                                   %
            if l > 0                                              %
                fprintf(self.FILE,'%s\n',...                      %
                    char(text(length(text))));                    %
            end                                                   %
        end                                                       %
                                                                  %
%% Write emty line to file ---------------------------------------%
        function writee(self)                                     %
            fprintf( self.FILE,' \n');                            %
        end                                                       %
                                                                  %
%% Close file ----------------------------------------------------%
        function finalize(self)                                   %
            fclose(self.FILE);                                    %
        end                                                       %
                                                                  %
    end                                                           %
                                                                  %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

