global files x y

% files = 'res\\lam-x%d.csv';
% x     = 2;
% y     = 9;
% format = 'k-';
% 
% main2();

steps = 1:3:16;

files = 'res\\turb0-x%d.csv';
x     = 2;
y     = 14;
format = 'k-*';

main2();

files = 'res\\turb3-x%d.csv';
x     = 2;
y     = 14;
format = 'k-s';

main2();
% 
% files = 'res\\turb4-x%d.csv';
% x     = 2;
% y     = 14;
% format = 'k-o';
% 
% main();
% 
% subplot(2,3,1);
% legend('lam','turb 0','turb 3','turb 4','Location','west');
% set(gcf,'color','w');