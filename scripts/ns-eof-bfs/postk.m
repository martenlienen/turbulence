%% Load CSV
f = csvread(file);
[h,w] =size(f);
Y = [];

%% Determin modes for each position on the path
for i=2:w
   
    x = f(1    ,i);
    y = f( max(2,end-la):end,i);
    
    Y(end+1,1) = x*1000;
    Y(end  ,2) = mean(y);
    Y(end  ,3) = var(y);
    Y(end  ,4) = skewness(y);
    Y(end  ,5) = kurtosis(y);
end

%% Modify output
i1 = find(Y(:,1)>mi,1,'first');
i2 = find(Y(:,1)<ma,1,'last' );
Y=Y(i1:i2,:);

Y(:,1) = fliplr(Y(:,1)')';
Y(:,1)=Y(:,1)-min(Y(:,1));

%% Plot
subplot(2,2,1)
plot(Y(:,1),Y(:,2),st,'LineWidth',2); grid on;hold on;
title('Mean')
xlabel('Distance from symmetry plane [mm]');
ylabel('Means [m/s]');

if stoch
    
    subplot(2,2,2)
    plot(Y(:,1),Y(:,3),st,'LineWidth',2); grid on;hold on;
    title('Variance')
    xlabel('Distance from symmetry plane [mm]');
    ylabel('Variance [m/s]');

    subplot(2,2,3)
    plot(Y(:,1),Y(:,4),st,'LineWidth',2); grid on;hold on;
    title('Skewness')
    xlabel('Distance from symmetry plane [mm]');
    ylabel('Skewness [-]');

    subplot(2,2,4)
    plot(Y(:,1),Y(:,5),st,'LineWidth',2); grid on;hold on;
    title('Kurtosis')
    xlabel('Distance from symmetry plane [mm]');
    ylabel('Kurtosis [-]');

end

set(gcf,'color','w');