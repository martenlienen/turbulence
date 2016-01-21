%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%          Function to visualize data in bode diagram           %%
%-----------------------------------------------------------------%
%                  File: BodeDataToBodeDiagram.m                  %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
function BodeDataToBodeDiagram( W, A, P )                         %
                                                                  %
    figure(10)                                                    %
                                                                  %
    %% Plot magnitude --------------------------------------------%
    subplot(2,1,1)                                                %
    semilogx(W,mag2db(A));                                        %
    hold on; grid on                                              %
    semilogx(W,mag2db(A),'x')                                     %
    xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');           %
                                                                  %
    %% Plot phase ------------------------------------------------%
    subplot(2,1,2)                                                %
    semilogx(W,P);                                                %
    hold on; grid on                                              %
    semilogx(W,P,'x');                                            %
    xlabel('Frequency [Hz]'); ylabel('Phase [deg]');              %
                                                                  %
    set(gcf,'color','w');                                         %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%