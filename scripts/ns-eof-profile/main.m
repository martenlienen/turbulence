
steps = 1:3:16;
counter = 1;
for s = steps
    
    
    subplot(2,3,counter);
    title(sprintf('x = %dm',s));
    sprintf(files,s)
    [X,Y] = getdata(sprintf(files,s),x,y);
    plot(X,Y,format); grid on;  hold on;
    
    xlabel('u [m/s]');
    ylabel('y [m]');
    
    counter = counter+1;
end

