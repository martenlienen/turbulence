clear all

global file Re l form colt lplus yplusb

%% Set up
l = 0.5;
yplusb = true;

%% Turbulence Model 0: no wall function

colt = reducedatacols('keturb-de');   colt.BY = colt.H;
file = 'data_2_13750_fine.csv';
Re   = 13750;                         lplus   = 1.333489e-03;
form = 'k';
tkemain();

axis([0 50 -0.03 0.04]);
xlabel('y+ [-]')
ylabel('TKE: gain and loss [...]');
title('TKE: gain and loss for Re=13750')
