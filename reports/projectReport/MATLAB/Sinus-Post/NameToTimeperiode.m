%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%         Function to extrect timeperiode from file name        %%
%-----------------------------------------------------------------%
%                    File:  NameToTimeperiode.m                   %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
function [ period ] = NameToTimeperiode( LOC )                    %
                                                                  %
    C = strsplit(LOC,'\\');                                       %
    C = C(length(C));                                             %
    LOCt = strrep(C(length(C)), 'res', '');                       %
    LOCt = strrep(LOCt, '.txt', '');                              %
    period = str2double(LOCt);                                    %
                                                                  %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

