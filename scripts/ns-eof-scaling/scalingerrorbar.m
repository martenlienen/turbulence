function [  ] = scalingerrorbar( xvalues, yvalues )

    xvaluestemp = [];
    yvaluestemp = [];
    errortemp = [];
    

    temp = [1; xvalues(1:end-1)<xvalues(2:end)]
    
    for i = 1:length(temp)
       if temp(i) == 1
            xvaluestemp(end+1) = xvalues(i);
            yvaluestemp(end+1) = yvalues(i);
            errortemp(end+1) = 0;
       else 
           a = (yvaluestemp(end)+yvalues(i))/2;
           e = abs(yvaluestemp(end)-yvalues(i));
            yvaluestemp(end) = a;
            errortemp(end) = e/2;
       end
        
    end
    
    errorbar(xvaluestemp,yvaluestemp,errortemp,'k');hold on; grid on;
    plot(xvaluestemp,yvaluestemp,'k*');


end

