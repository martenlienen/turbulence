%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                  Function to visualize errors                 %%
%-----------------------------------------------------------------%
%                     File:  BodeDataToError.m                    %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
function  BodeDataToError( W, A, P )                              %
                                                                  %
    figure(40)                                                    %
                                                                  %
    %% Error for the input signal --------------------------------%
    subplot(2,1,1)                                                %
    semilogx(W,A*100);                                            %
    hold on; grid on;                                             %
    semilogx(W,A*100,'x')                                         %
    xlabel('Frequency [Hz]');   ylabel('Error IN [%]')            %
                                                                  %
    %% Error for the output signal -------------------------------%
    subplot(2,1,2)                                                %
    semilogx(W,P*100);                                            %
    hold on; grid on                                              %
    semilogx(W,P*100,'x');                                        %
    xlabel('Frequency [Hz]');   ylabel('Error OUT [%]')           %
                                                                  %
    set(gcf,'color','w');                                         %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%