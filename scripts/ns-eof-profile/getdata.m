function [ U, Y ] = getdata( file, c1,c2 )

M = csvread(file,1,0);
U = M(1:end,c1);
Y = M(1:end,c2);

end

