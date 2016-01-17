%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                               %%
%-----------------------------------------------------------------%
%                           File: gf2.m                           %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
%% Initialize ----------------------------------------------------%
A = RES(:,5);                                                     %
P = RES(:,6);                                                     %
F = RES(:,4);                                                     %
                                                                  %
%% Get gain function ---------------------------------------------%
% Here you have to run System Identifcation Tool using...         %
%      A as magnitude                                             %
%      P as periode                                               %
%      F as frequency                                             %
                                                                  %
%% Plot the gain function in the bode diagaram with orig. signal -%
[A1,P1]=bode(F,tf1);                                              %
BodeDataToBodeDiagram(F/2/pi,A,P)                                 %
BodeDataToBodeDiagram...                                          %
           (F/2/pi,squeeze(A1(1,1,:)),squeeze(P1(1,1,:)))         %
legend('experiment','experiment - v','model','model - v')         %
                                                                  %
%% Plot diagrams for gain function -------------------------------%
figure(50)                                                        %
subplot(2,2,1)                                                    %
bode(tf1);grid on;                                                %
subplot(2,2,2)                                                    %
nyquist(tf1);grid on;                                             %
subplot(2,2,3)                                                    %
rlocus(tf1);grid on;                                              %
subplot(2,2,4)                                                    %
nichols(tf1);grid on;                                             %
set(gcf,'color','w');                                             %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%