global files x y

files = 'res\\lam-x%d.csv';
x     = 2;
y     = 9;
format = 'k-';

main();

files = 'res\\turb0-x%d.csv';
x     = 2;
y     = 14;
format = 'k-*';

main();

files = 'res\\turb3-x%d.csv';
x     = 2;
y     = 14;
format = 'k-s';

main();

files = 'res\\turb4-x%d.csv';
x     = 2;
y     = 14;
format = 'k-o';

main();

subplot(2,3,1);
legend('lam','turb 0','turb 3','turb 4','Location','west');
set(gcf,'color','w');