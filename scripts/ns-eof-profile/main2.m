% files = 'res\\turb0-x%d.csv';
% x     = 2;
% y     = 14;
% format = 'k-*';



for i = steps
    xtemp = i;
    
    [X,Y] = getdata(sprintf(files,i),x,y);
    plot(xtemp+X*scale,Y,format);
    grid on; hold on;
    plot([xtemp xtemp],[min(Y) max(Y)] ,'k');
    

end

set(gcf,'color','w');
xlabel('x [m]');
ylabel('y [m]');