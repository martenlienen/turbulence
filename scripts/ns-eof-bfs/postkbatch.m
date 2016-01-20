global file
	
%% setup
% st    = style
% stoch = should variace, skewness, kurtosisi be calculated
% la    = only consider the last la timesteps
% mi    = lower x-positon
% ma    = higher x-position
	
st = 'k'; stoch = true;
la = 500;
mi = 5; ma = 105;

file = 'karman350.csv';

%% Run
postk();
