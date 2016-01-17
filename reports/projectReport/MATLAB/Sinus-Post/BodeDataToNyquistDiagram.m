%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%        Function to visualize data in nyquist diagram          %%
%-----------------------------------------------------------------%
%                 File: BodeDataToNyquistDiagram.m                %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
function BodeDataToNyquistDiagram( W, A, P )                      %
                                                                  %
    % Amplitude and phase to pointer                              %
    p=P*pi/180;     a=A.*(cos(p)+1i*sin(p));                      %
                                                                  %
    %% Plot Nyquist (cartesian coordinate system) ----------------%
    figure(20)                                                    %
    plot(a); hold on; plot(a,'x'); grid on;                       %
    xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');           %
    set(gcf,'color','w');                                         %
                                                                  %
    %% Plot Nyquist (cylindical coordinate system) ---------------%
    figure(30)                                                    %
    temp=0:0.01:2*pi;                                             %
    polar(0:0.01:2*pi,temp./temp);hold on;                        %
    polar(p,A);  polar(p,A,'x');                                  %
                                                                  %
    set(gcf,'color','w');                                         %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

