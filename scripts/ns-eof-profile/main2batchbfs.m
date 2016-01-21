global files x y


steps = 0.1:0.1:1.1;

files = 'bfs\\aturb%.1f.csv';
x     = 2;
y     = 15;
format = 'b';
scale = 0.05;
main2();

files = 'bfs\\keturb%.1f.csv';
x     = 2;
y     = 23;
format = 'r';
% scale = 0.1;
main2();

axis([0 1.25 0 0.2]);