%% Load relevant data at the centre of each cell
%  Y  - Y-position
%  U  - velocity in x-direction
%  NU - eddy viscosity
%  US - velocity fluctuation

[U,  Y] = reducedata(file,colt.BY,colt.ARCL,colt.U1, div);
[NU, Y] = reducedata(file,colt.BY,colt.ARCL,colt.NUT,div);
[US, Y] = reducedata(file,colt.BY,colt.ARCL,colt.u  ,div);


%% Calculate help variables
%  S12 - shear stress

S12 = ((U(3:end)-U(1:end-2))./(Y(3:end)-Y(1:end-2)));

%% Recalculate variables (algebraic model) and compare to measured values
%  RNU - eddy viscosity
%  RUS - velocity fluctuation

% RNU = (0.41*Y(2:end-1)).^(2).*sqrt(2*S12.*S12);
% RUS = (0.41*Y(2:end-1)).*sqrt(4*S12.*S12);
% 
% figure(1000)
% plot(Y(2:end-1),RNU);hold on;grid on;
% plot(Y,NU);
% 
% figure(1100)
% plot(Y(2:end-1),RUS);hold on;grid on;
% plot(Y,US);

%% Calculate shear stresses
% X       - wall distance
% taul(X) - laminar shear stress
% taut(X) - reynolds shear stress
% tauw    - wall shear stress (taul+taut at the posion X=0)

X = Y(2:end-1);
taul = S12/Re;
taut = S12.*NU(2:end-1);

p = polyfit(X,taul+taut,1);

if(tauwb)
tauw = polyval(p,0);
end
%% Calculate ...
% utau  - ...
% lplus - ...

    utau = sqrt(tauw);


lplus = 1/utau/Re;

disp(sprintf('Re = %4d: l+ = %.6d tauw = %.6d utau = %.6d',Re,lplus,tauw,utau))

%% Plot shear stresses

figure(2000)

%left: laminar shear stress over wall distance
subplot(1,3,1)
plot([0; X; 0.5]/l,[1; taul/tauw; 0],form,'LineWidth',2);
grid on;hold on;
% plot([0; X; 0.5]/l,[1; (taut+taul)/tauw; 0],form,'LineWidth',1);

ylabel('\tau_v/\tau_w [-]');
xlabel('1-y/h [-]');

% center: reynolds shear stress over wall distance
subplot(1,3,2)
plot([0; X; 0.5]/l,[0; taut/tauw; 0],form,'LineWidth',2);
grid on;hold on;
% plot([0; X; 0.5]/l,[1; (taut+taul)/tauw; 0],form,'LineWidth',1);

ylabel('\tau_r/\tau_w [-]');
xlabel('1-y/h [-]');

% right: laminar & reynolds shear stress over y+
subplot(1,3,3)
plot([0; X; 0.5]/lplus,[0; taut./(taut+taul); 0],form,'LineWidth',2);grid on;hold on;
plot([0; X; 0.5]/lplus,[1; taul./(taut+taul); 0],form,'LineWidth',2);

ylabel('\tau_v/\tau_w [-] and \tau_t/\tau_w ');
xlabel('y/y+ [-]');

xlim([0 50]);

%% Plot w+ over y+
figure(3000)
semilogx([0; Y]/lplus,[0; U]/utau,form,'LineWidth',2);grid on; hold on;

